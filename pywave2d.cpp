#include <Python.h>
#include <stdio.h>
#include <numpy/arrayobject.h>

#include "pywave2d.hpp"
#include "pywave2dexception.hpp"

static const double DEFAULT_PB = 0.05;
static const int DEFAULT_SEED = 983;

static const int WAVELET_SUBTILE = 1;

const double pyWave2D::INVALID_PROPOSAL = -1.0e30;

static const pyWave2D::wavelet2d_inverse_t WAVELET_INVERSE[] = {
  haar_lift_inverse2d_haar,
  daub4_dwt_inverse2d_daub4,
  daub6_dwt_inverse2d_daub6,
  daub8_dwt_inverse2d_daub8,
  cdf97_lift_inverse2d_cdf97,
  cdf97_lift_periodic_inverse2d_cdf97
};

static int coordtoindex(void *user, int i, int j, int k, int depth);
static int indextocoord(void *user, int index, int *i, int *j, int *k, int *depth);


pyWave2D::pyWave2D(int degreew,
		   int degreeh,
		   int wavelet,
		   int _kmax) :
  wavelet2d_inverse(WAVELET_INVERSE[0]),
  Pb(0.05),
  width(-1),
  height(-1),
  size(-1),
  ncoeff(-1),
  wt(NULL),
  proposal(NULL),
  coefficient_histogram(NULL),
  kmax(_kmax),
  hnk(NULL),
  rng(NULL),
  seed(DEFAULT_SEED),
  model(NULL),
  workspace(NULL),
  last_proposal(-1),
  last_level(-1),
  pbirth(NULL),
  abirth(NULL),
  pdeath(NULL),
  adeath(NULL),
  pvalue(NULL),
  avalue(NULL),
  fb(-1),
  fd(-1),
  pb(-1),
  pd(-1),
  indices(new int[_kmax]),
  values(new double[_kmax])
{
  if (degreew < 1 || degreew > 16) {
    throw PYWAVE2DEXCEPTION("Invalid degree w %d", degreew);
  }

  if (degreeh < 1 || degreeh > 16) {
    throw PYWAVE2DEXCEPTION("Invalid degree h %d", degreeh);
  }

  if (wavelet >= 0 && wavelet < (int)(sizeof(WAVELET_INVERSE)/sizeof(wavelet2d_inverse_t))) {
    wavelet2d_inverse = WAVELET_INVERSE[wavelet];
  } else {
    throw PYWAVE2DEXCEPTION("Invalid wavelet parameter");
  }

  wt = wavetree2d_sub_create(degreew, degreeh, 0.0);
  if (wt == NULL) {
    throw PYWAVE2DEXCEPTION("Failed to create wavetree");
  }

  width = wavetree2d_sub_get_width(wt);
  height = wavetree2d_sub_get_height(wt);
  size = wavetree2d_sub_get_size(wt);
  ncoeff = wavetree2d_sub_get_ncoeff(wt);
  maxdepth = wavetree2d_sub_maxdepth(wt);
  
  hnk = hnk_cartesian_nonsquare_2D_create_sub(degreew,
					      degreeh,
					      kmax);
  if (hnk == NULL) {
    throw PYWAVE2DEXCEPTION("Failed to create hnk");
  }

  rng = gsl_rng_alloc(gsl_rng_taus);
  if (rng == NULL) {
    throw PYWAVE2DEXCEPTION("Failed to allocate random number generator\n");
  }
  
  gsl_rng_set(rng, seed);

  model = new double[size];
  int workspace_size = width;
  if (height > width) {
    workspace_size = height;
  }
  workspace = new double[workspace_size];

  coefficient_histogram = coefficient_histogram_create(ncoeff,
						       100,
						       -1.0,
						       1.0,
						       coordtoindex,
						       indextocoord,
						       wt);
  if (coefficient_histogram == NULL) {
    throw PYWAVE2DEXCEPTION("Failed to create coefficient histogram");
  }

  Pb = DEFAULT_PB;
  last_proposal = WT_PERTURB_INVALID;
  last_level = -1;

  pbirth = new int[maxdepth + 1];
  abirth = new int[maxdepth + 1];
  pdeath = new int[maxdepth + 1];
  adeath = new int[maxdepth + 1];
  pvalue = new int[maxdepth + 1];
  avalue = new int[maxdepth + 1];

  for (int i = 0; i <= maxdepth; i ++) {
    pbirth[i] = 0;
    abirth[i] = 0;
    pdeath[i] = 0;
    adeath[i] = 0;
    pvalue[i] = 0;
    avalue[i] = 0;
  }
		   
}

pyWave2D::~pyWave2D()
{
}

void
pyWave2D::set_log(const char *filename)
{
  if (filename != NULL) {
    slog_set_output_file(filename, SLOG_FLAGS_CLEAR);
  }
}

void
pyWave2D::set_seed(int _seed)
{
  seed = _seed;
  gsl_rng_set(rng, seed);
}

bool
pyWave2D::load_priorproposal(const char *filename)
{
  proposal = wavetree_pp_load(filename, seed, coefficient_histogram);
  if (proposal == NULL) {
    ERROR("fwave2d_create: failed to load/create proposal: \"%s\"\n", filename);
    return false;
  }

  for (int i = 0; i < ncoeff; i ++) {
    int depth = wavetree2d_sub_depthofindex(wt, i);

    int ii;
    int ij;

    double vmin;
    double vmax;
      
    if (wavetree2d_sub_2dindices(wt, i, &ii, &ij) < 0) {
      ERROR("fwave2d_create: failed to determine 2d indices\n");
      return false;
    }

    if (wavetree_pp_prior_range2d(proposal,
                                  ii,
                                  ij,
				  depth,
                                  maxdepth,
                                  0.0,
                                  &vmin,
                                  &vmax) < 0) {
      ERROR("fwave2d_create: failed to get coefficient range\n");
      return false;
    }

    if (coefficient_histogram_set_range(coefficient_histogram,
                                        i,
                                        vmin,
                                        vmax) < 0) {
      ERROR("fwave2d_create: failed to set range\n");
      return false;
    }
  }

  return true;
}
  
  
void
pyWave2D::set_pb(double pb)
{
  if (pb < 0.0 || pb > 0.5) {
    throw PYWAVE2DEXCEPTION("Pb out of range: %f", pb);
  }

  Pb = pb;
}

void
pyWave2D::initialise(double dc)
{
  if (wavetree2d_sub_initialize(wt, dc) < 0) {
    ERROR("Failed to initialise wave tree");
  }
}

double
pyWave2D::propose()
{
  if (last_proposal != WT_PERTURB_INVALID) {
    throw PYWAVE2DEXCEPTION("Perturbation exists");
  }
  
  if (proposal == NULL) {
    throw PYWAVE2DEXCEPTION("No prior/proposal loaded");
  }

  double u = gsl_rng_uniform(rng);
  int r = -1;

  double ppratio = 0.0;
  
  if (u < Pb) {

    pb ++;
   
    /*
     * Birth
     */
    last_proposal = WT_PERTURB_BIRTH;
    r = birth(kmax,
	      maxdepth,
	      wt,
	      proposal,
	      hnk,
	      rng,
	      &ppratio,
	      &last_level);
    if (last_level >= 0) {
      pbirth[last_level] ++;
    } else {
      fb ++;
    }
    if (r < 0) {
      ERROR("Birth error");
      last_proposal = WT_PERTURB_INVALID;
      return INVALID_PROPOSAL;
    } else if (r == 0) {
      last_proposal = WT_PERTURB_INVALID;
      return INVALID_PROPOSAL;
    }

  } else if (u < (2.0 * Pb)) {

    pd ++;
    
    /* 
     * Death
     */
    last_proposal = WT_PERTURB_DEATH;
    r = death(kmax,
	      maxdepth,
	      wt,
	      proposal,
	      hnk,
	      rng,
	      &ppratio,
	      &last_level);
    if (last_level >= 0) {
      pdeath[last_level] ++;
    } else {
      fd ++;
    }
    if (r < 0) {
      ERROR("Death error");
      last_proposal = WT_PERTURB_INVALID;
      return INVALID_PROPOSAL;
    } else if (r == 0) {
      last_proposal = WT_PERTURB_INVALID;
      return INVALID_PROPOSAL;
    }

  } else {
    
    /*
     * Value
     */
    last_proposal = WT_PERTURB_VALUE;
    r = value(kmax,
	      maxdepth,
	      wt,
	      proposal,
	      hnk,
	      rng,
	      &ppratio,
	      &last_level);
    if (last_level >= 0) {
      pvalue[last_level] ++;
    }
    
    if (r < 0) {
      ERROR("Value error");
      last_proposal = WT_PERTURB_INVALID;
      return INVALID_PROPOSAL;
    } else if (r == 0) {
      last_proposal = WT_PERTURB_INVALID;
      return INVALID_PROPOSAL;
    }
  }

  return ppratio;
}
  
void
pyWave2D::accept()
{
  if (last_proposal == WT_PERTURB_INVALID) {
    ERROR("No perturbation");
  }
  
  if (wavetree2d_sub_commit(wt) < 0) {
    ERROR("fwave2d_accept: failed to commit pertubation %d\n", last_proposal);
  }

  if (last_level < 0 || last_level > maxdepth) {
    ERROR("fwave2d_accept: invalid level for statistics %d\n", last_level);
  }
  
  switch (last_proposal) {
  case WT_PERTURB_BIRTH:
    abirth[last_level] ++;
    break;

  case WT_PERTURB_DEATH:
    adeath[last_level] ++;
    break;

  case WT_PERTURB_VALUE:
    avalue[last_level] ++;
    break;

  default:
    ERROR("fwave2d_accepth: unknown perturbation in statistics\n");
  }
  
  last_proposal = WT_PERTURB_INVALID;
}


void
pyWave2D::reject()
{
  if (last_proposal == WT_PERTURB_INVALID) {
    ERROR("No perturbation");
  }

  if (wavetree2d_sub_undo(wt) < 0) {
    ERROR("fwave2d_accept: failed to undo pertubation\n");
  }

  last_proposal = WT_PERTURB_INVALID;
}

int
pyWave2D::levels() const
{
  return maxdepth + 1;
}
  
void
pyWave2D::pb_counts(int *n, int **counts) const
{
  *n = maxdepth + 1;
  *counts = pbirth;
}

void
pyWave2D::ab_counts(int *n, int **counts) const
{
  *n = maxdepth + 1;
  *counts = abirth;
}

void
pyWave2D::pd_counts(int *n, int **counts) const
{
  *n = maxdepth + 1;
  *counts = pdeath;
}

void
pyWave2D::ad_counts(int *n, int **counts) const
{
  *n = maxdepth + 1;
  *counts = adeath;
}

void
pyWave2D::pv_counts(int *n, int **counts) const
{
  *n = maxdepth + 1;
  *counts = pvalue;
}

void
pyWave2D::av_counts(int *n, int **counts) const
{
  *n = maxdepth + 1;
  *counts = avalue;
}

int
pyWave2D::get_width() const
{
  return width;
}

int
pyWave2D::get_height() const
{
  return height;
}

void
pyWave2D::get_image(int *w, int *h, double **_image)
{
  memset(model, 0, sizeof(double) * size);
  if (wavetree2d_sub_map_to_array(wt, model, size) < 0) {
    ERROR("get_image: failed to map model to coefficient image\n");
    return;
  }

  if (wavelet2d_inverse(model,
			width,
			height,
			width,
			workspace,
			WAVELET_SUBTILE) < 0) {
    ERROR("get_image: failed to do inverse transform\n");
    return;
  }

  *w = width;
  *h = height;
  *_image = model;
}

void
pyWave2D::model_retrieve()
{
  if (wavetree2d_sub_get_model(wt, kmax, indices, values, &k) < 0) {
    throw PYWAVE2DEXCEPTION("Failed to get model");
  }
}

void
pyWave2D::get_model_indices(int *n, int **_indices) const
{
  *n = k;
  *_indices = indices;
}

void
pyWave2D::get_model_values(int *nv, double **_values) const
{
  *nv = k;
  *_values = values;
}

void
pyWave2D::set_model_indices(int n, const int *_indices)
{
  if (n > kmax) {
    throw PYWAVE2DEXCEPTION("Too many indices %d > %d", n, kmax);
  }

  k = n;
  for (int i = 0; i < k; i ++) {
    indices[i] = _indices[i];
  }
}

void
pyWave2D::set_model_values(int nv, const double *_values)
{
  if (nv != k) {
    throw PYWAVE2DEXCEPTION("Mismatch between indices and values: %d != %d", nv, k);
  }

  for (int i = 0; i < k; i ++) {
    values[i] = _values[i];
  }
}

void
pyWave2D::model_commit()
{
  if (wavetree2d_sub_set_model(wt, indices, values, k) < 0) {
    throw PYWAVE2DEXCEPTION("Failed to set model");
  }
}

int pyWave2D::get_k() const
{
  return wavetree2d_sub_coeff_count(wt);
}

double
pyWave2D::get_prior_scale() const
{
  double v;

  wavetree_pp_setscale(proposal, 0.0, &v);

  return v;
}

void
pyWave2D::set_prior_scale(double scale)
{
  wavetree_pp_setscale(proposal, scale, NULL);
} 

double
pyWave2D::propose_prior_scale(double new_scale)
{
  double current_log_prior = wavetree2d_sub_logpriorprobability(wt, proposal);
  double old_scale;

  wavetree_pp_setscale(proposal, new_scale, &old_scale);

  double proposed_log_prior = wavetree2d_sub_logpriorprobability(wt, proposal);

  wavetree_pp_setscale(proposal, old_scale, NULL);

  return proposed_log_prior - current_log_prior;
}

double
pyWave2D::get_logprior()
{
  return wavetree2d_sub_logpriorprobability(wt, proposal);
}

int
pyWave2D::death(int kmax,
		int maxdepth,
		wavetree2d_sub_t *wt,
		wavetree_pp_t *proposal,
		hnk_t *hnk,
		gsl_rng *rng,
		double *priorpropratio,
		int *level)
{
  int k;

  double ratio;

  int death_depth;
  int death_idx;
  double death_value;
  double death_parent_coeff;
  
  double choose_prob;
  double reverse_prob;
  double death_prob;
  double prior_prob;
  
  int ii;
  int ij;

  k = wavetree2d_sub_coeff_count(wt);
  *level = -1;
  
  if (k > 1) {

    if (hnk_get_kplus1_ratio(hnk,
                             maxdepth,
                             k - 1,
                             &ratio) < 0) {
      ERROR("failed to get ratio fo death");
      return -1;
    }
    
    if (wavetree2d_sub_choose_death_global(wt,
                                           gsl_rng_uniform(rng),
                                           maxdepth,
                                          &death_depth,
                                           &death_idx,
                                           &choose_prob) < 0) {
      ERROR("failed to choose death index");
      return -1;
    }
    *level = death_depth;
    
    if (wavetree2d_sub_propose_death(wt,
                                     death_idx,
                                     death_depth,
                                     &death_value) < 0) {
      ERROR("failed to propose death");
      return -1;
    }

    if (wavetree2d_sub_2dindices(wt, death_idx, &ii, &ij) < 0) {
      ERROR("failed to compute 2d indices");
      return -1;
    }

    if (wavetree2d_sub_get_coeff(wt,
                                 wavetree2d_sub_parent_index(wt, death_idx),
                                 &death_parent_coeff) < 0) {
      ERROR("failed to get parent coefficient");
      return -1;
    }

    if (wavetree_pp_death2d(proposal,
                            ii,
                            ij,
                            death_depth,
                            maxdepth,
                            death_parent_coeff,
                            death_value,
                            &death_prob) < 0) {
      ERROR("failed to get death probability");
      return -1;
    }

    if (wavetree2d_sub_reverse_death_global(wt,
                                            maxdepth,
                                            death_depth,
                                            death_idx,
                                            &reverse_prob) < 0) {
      ERROR("failed to do reverse death");
      return -1;
    }

    prior_prob = wavetree_pp_prior_probability2d(proposal,
                                                 ii, ij,
                                                 death_depth,
                                                 maxdepth,
                                                 death_parent_coeff,
                                                 death_value);

    /* printf("%10.6f %10.6f %10.6f %10.6f %10.6f\n", */
    /*     reverse_prob, */
    /*     choose_prob, */
    /*     death_prob, */
    /*     ratio, */
    /*     prior_prob); */

    *priorpropratio =
      log(reverse_prob) -
      log(choose_prob) +
      log(death_prob) +
      log(ratio) -
      log(prior_prob);

    return 1;
  }

  return 0;
}


int
pyWave2D::birth(int kmax,
		int maxdepth,
		wavetree2d_sub_t *wt,
		wavetree_pp_t *proposal,
		hnk_t *hnk,
		gsl_rng *rng,
		double *priorpropratio,
		int *level)
{
  int k;
  double ratio;

  int birth_depth;
  int birth_idx;
  double birth_value;
  double birth_parent_coeff;

  double choose_prob;
  double reverse_prob;
  double birth_prob;
  double prior_prob;
  
  int birth_valid;

  int ii;
  int ij;

  *level = -1;

  k = wavetree2d_sub_coeff_count(wt);
  if (k < 0) {
    ERROR("failed to get no. coefficients.");
    return -1;
  }
  
  if (k == kmax) {
    return 0;
  }

  if (hnk_get_kplus1_ratio(hnk, maxdepth, k, &ratio) < 0) {
    ERROR("failed to compute hnk ratio");
    return -1;
  }

  if (ratio <= 0.0) {
    ERROR("invalid ratio %f for k %d (kmax %d, maxdepth %d)", ratio, k, kmax, maxdepth);
    return -1;
  }

  ratio = 1.0/ratio;

  if (wavetree2d_sub_choose_birth_global(wt,
                                         gsl_rng_uniform(rng),
                                         maxdepth,
                                         &birth_depth,
                                         &birth_idx,
                                         &choose_prob) < 0) {
    ERROR("failed to choose birth");
    return -1;
  }
  *level = birth_depth;

  birth_valid = 0;

  if (wavetree2d_sub_2dindices(wt, birth_idx, &ii, &ij) < 0) {
    ERROR("failed to compute 2d indices");
    return -1;
  }

  if (wavetree2d_sub_get_coeff(wt,
                               wavetree2d_sub_parent_index(wt, birth_idx),
                               &birth_parent_coeff) < 0) {
    ERROR("failed to get parent coefficient");
    return -1;
  }

  if (wavetree_pp_birth2d(proposal,
                          ii,
                          ij,
                          birth_depth,
                          maxdepth,
                          birth_parent_coeff,
                          &birth_value,
                          &birth_prob,
                          &birth_valid) < 0) {
    ERROR("failed to do proposal");
    return -1;
  }

  if (birth_valid) {

    if (wavetree2d_sub_propose_birth(wt,
                                     birth_idx,
                                     birth_depth,
                                     birth_value) < 0) {
      ERROR("failed to propose birth");
      return -1;
    }

    if (wavetree2d_sub_reverse_birth_global(wt,
                                            maxdepth,
                                            birth_depth,
                                            birth_idx,
                                            &reverse_prob) < 0) {
      ERROR("failed to reverse birth");
      return -1;
    }

    prior_prob = wavetree_pp_prior_probability2d(proposal,
                                                 ii,
                                                 ij,
                                                 birth_depth,
                                                 maxdepth,
                                                 birth_parent_coeff,
                                                 birth_value);


    /* printf("%10.6f %10.6f %10.6f %10.6f %10.6f\n", */
    /*     reverse_prob, */
    /*     choose_prob, */
    /*     birth_prob, */
    /*     ratio, */
    /*     prior_prob); */
    
    *priorpropratio =
      log(reverse_prob) -
      log(choose_prob) -
      log(birth_prob) +
      log(ratio) +
      log(prior_prob);
    
    return 1;
  }

  return 0;
}

int
pyWave2D::value(int kmax,
		int maxdepth,
		wavetree2d_sub_t *wt,
		wavetree_pp_t *proposal,
		hnk_t *hnk,
		gsl_rng *rng,
		double *priorpropratio,
		int *level)
{
  int value_depth;
  int value_idx;
  int parent_idx;
  double value;
  double value_parent_coeff;
  
  double choose_prob;

  int ii;
  int ij;

  double value_prior_ratio;
  int prior_errors;
  double temperature = 1.0;

  if (wavetree2d_sub_choose_value_global(wt,
                                         gsl_rng_uniform(rng),
                                         maxdepth,
                                         &value_depth,
                                         &value_idx,
                                         &choose_prob) < 0) {
    ERROR("failed to choose value");
    return -1;
  }
  
  *level = value_depth;

  if (wavetree2d_sub_get_coeff(wt, value_idx, &value) < 0) {
    ERROR("failed to get coefficient value");
    return -1;
  }

  if (wavetree2d_sub_2dindices(wt, value_idx, &ii, &ij) < 0) {
    ERROR("failed to get 2d indices for coefficient");
    return -1;
  }

  if (wavetree_pp_value_init(proposal) < 0) {
    ERROR("failed to initialise proposal");
    return -1;
  }

  if (value_idx > 0) {
    parent_idx = wavetree2d_sub_parent_index(wt, value_idx);
    if (parent_idx < 0 ) {
      ERROR("failed to get parent index.");
      fprintf(stderr, "no parent for %d\n", value_idx);
      return -1;
    }
    if (wavetree2d_sub_get_coeff(wt,
                                 parent_idx,
                                 &value_parent_coeff) < 0) {
      ERROR("failed to get parent coefficient");
      return -1;
    }
  } else {
    parent_idx = -1;
    value_parent_coeff = 0.0;
  }

  value_prior_ratio = 1.0;
  if (wavetree_pp_propose_value2d(proposal,
                                  ii,
                                  ij,
                                  value_depth,
                                  maxdepth,
                                  value_parent_coeff,
				  temperature,
                                  &value,
                                  &value_prior_ratio) < 0) {
    ERROR("failed to perturb value");
    return -1;
  }

  prior_errors = wavetree_pp_value_error_count(proposal);
  if (prior_errors < 0) {
    ERROR("failed to check errors");
    return -1;
  }

  if (prior_errors == 0) {
    if (wavetree2d_sub_propose_value(wt, value_idx, value_depth, value) < 0) {
      ERROR("failed to propose new value");
      return -1;
    }

    *priorpropratio = log(value_prior_ratio);

    return 1;
  }

  return 0;
}

static int coordtoindex(void *user, int i, int j, int k, int depth)
{
  wavetree2d_sub_t *wt = (wavetree2d_sub_t*)user;

  return wavetree2d_sub_from_2dindices(wt, i, j);
}

static int indextocoord(void *user, int index, int *i, int *j, int *k, int *depth)
{
  wavetree2d_sub_t *wt = (wavetree2d_sub_t*)user;

  if (wavetree2d_sub_2dindices(wt, index, i, j) < 0) {
    return -1;
  }

  *k = 0;
  *depth = wavetree2d_sub_depthofindex(wt, index);
  return 0;
}
