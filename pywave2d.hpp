//
//    pyWave2D : A python 2 extension for inversion of 2D datasets using trans-dimensional
//    trees. For reference see
//
//      R Hawkins and M Sambridge, "Geophysical imaging using trans-dimensional trees",
//      Geophysical Journal International, 2015, 203:2, 972 - 1000,
//      https://doi.org/10.1093/gji/ggv326
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//
#pragma once
#ifndef pywave2d_hpp
#define pywave2d_hpp

#include <gmp.h>
#include <gsl/gsl_rng.h>

extern "C" {
  
#include "slog.h"

#include "wavetree2d_sub.h"
#include "hnk_cartesian_nonsquare.h"

#include "cdf97_lift.h"
#include "cdf97_lift_periodic.h"
#include "haar_lift.h"
#include "daub4_dwt.h"
#include "daub6_dwt.h"
#include "daub8_dwt.h"
};

enum {
  WAVELET_HAAR = 0,
  WAVELET_DAUB4 = 1,
  WAVELET_DAUB6 = 2,
  WAVELET_DAUB8 = 3,
  WAVELET_CDF97 = 4,
  WAVELET_CDF97_PERIODIC = 5
};

class pyWave2D {
public:

  static const double INVALID_PROPOSAL;
  

  typedef int (*wavelet2d_inverse_t)(double *s,
				     int width,
				     int height,
				     int stride,
				     double *work,
				     int subtile);

  pyWave2D(int degreew,
	   int degreeh,
	   int wavelet,
	   int kmax);
  ~pyWave2D();

  void set_log(const char *filename);
  
  void set_seed(int seed);

  bool load_priorproposal(const char *filename);
  
  void set_pb(double pb);

  void initialise(double dc);
  
  double propose();
  
  void accept();

  void reject();

  int levels() const;
  
  void pb_counts(int *n, int **counts) const;
  void ab_counts(int *n, int **counts) const;
  void pd_counts(int *n, int **counts) const;
  void ad_counts(int *n, int **counts) const;
  void pv_counts(int *n, int **counts) const;
  void av_counts(int *n, int **counts) const;

  void get_image(int *w, int *h, double **image);

  int get_width() const;

  int get_height() const;

  int get_k() const;

  double get_prior_scale() const;

  void set_prior_scale(double scale);

  double propose_prior_scale(double new_scale);

  double get_logprior();
  
  void model_retrieve();
  void get_model_indices(int *n, int **indices) const;
  void get_model_values(int *nv, double **values) const;

  void set_model_indices(int n, const int *indices);
  void set_model_values(int nv, const double *values);
  void model_commit();

private:

  int death(int kmax,
	    int maxdepth,
	    wavetree2d_sub_t *wt,
	    wavetree_pp_t *proposal,
	    hnk_t *hnk,
	    gsl_rng *rng,
	    double *priorpropratio,
	    int *level);

  int birth(int kmax,
	    int maxdepth,
	    wavetree2d_sub_t *wt,
	    wavetree_pp_t *proposal,
	    hnk_t *hnk,
	    gsl_rng *rng,
	    double *priorpropratio,
	    int *level);

  int value(int kmax,
	    int maxdepth,
	    wavetree2d_sub_t *wt,
	    wavetree_pp_t *proposal,
	    hnk_t *hnk,
	    gsl_rng *rng,
	    double *priorpropratio,
	    int *level);


  wavelet2d_inverse_t wavelet2d_inverse;

  double Pb;

  int width;
  int height;
  int size;
  int ncoeff;
  wavetree2d_sub_t *wt;

  wavetree_pp_t *proposal;
  coefficient_histogram_t *coefficient_histogram;

  int kmax;
  int maxdepth;

  hnk_t *hnk;

  gsl_rng *rng;
  int seed;

  double *model;
  double *workspace;

  int last_proposal;
  int last_level;

  int *pbirth;
  int *abirth;

  int *pdeath;
  int *adeath;

  int *pvalue;
  int *avalue;

  int fb;
  int fd;
  
  int pb;
  int pd;

  int k;
  int *indices;
  double *values;

  
};

#endif // pywave2d_hpp
