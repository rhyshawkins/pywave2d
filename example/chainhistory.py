

import pickle

class ChainHistory:

    DELTA_INITIAL = 0
    DELTA_ADD = 1
    DELTA_REMOVE = 2
    DELTA_CHANGE = 3

    def __init__(self):

        self.ch = []
        self.model = None

    def __len__(self):
        return len(self.ch)

    def __iter__(self):

        (like, norm, errscale, priorscale), (code, current_indices, current_values) = self.ch[0]
        if code != ChainHistory.DELTA_INITIAL:
            raise Exception('Unexpected first entry in history %d' % self.ch[0][0])

        yield ((like, norm, errscale, priorscale), current_indices, current_values)
        
        for (like, norm, errscale, priorscale), h in self.ch[1:]:

            for (code, idx, val) in h:

                if (code == ChainHistory.DELTA_INITIAL):
                    current_indices = idx
                    current_values = val

                elif (code == ChainHistory.DELTA_ADD):
                    current_indices, current_values = self.__add(current_indices, current_values, idx, val)

                elif (code == ChainHistory.DELTA_REMOVE):
                    current_indices, current_values = self.__remove(current_indices, current_values, idx, val)

                elif (code == ChainHistory.DELTA_CHANGE):
                    current_indices, current_values = self.__change(current_indices, current_values, idx, val)

                else:
                    raise Exception('Unknown code %d' % code)

            yield ((like, norm, errscale, priorscale), current_indices, current_values)
                    
    def add_initial(self, like, norm, errscale, priorscale, indices, values):

        self.ch = [((like, norm, errscale, priorscale), (ChainHistory.DELTA_INITIAL, list(indices), list(values)))]

    def add_delta(self, like, norm, errscale, priorscale, old_indices, old_values, new_indices, new_values):

        i0 = list(old_indices)
        v0 = list(old_values)
        iv0 = zip(i0, v0)
        iv0.sort()
        i0, v0 = zip(*iv0)
        
        i1 = list(new_indices)
        v1 = list(new_values)
        iv1 = zip(i1, v1)
        iv1.sort()
        i1, v1 = zip(*iv1)
        
        self.ch.append(((like, norm, errscale, priorscale), (self.__compute_difference(i0, v0, i1, v1))))

    def save(self, filename):

        f = open(filename, 'w')
        pickle.dump(self.ch, f)
        f.close()

    def load(self, filename):

        f = open(filename, 'r')
        self.ch = pickle.load(f)
        self.model = None
        f.close()

    def dump(self, filename):

        f = open(filename, 'w')
        
        for h in self.ch:

            if type(h) == list:
                for (code, idx, val) in h:
                    
                    f.write('[%d %d %f] ' % (code, idx, val))

                f.write('\n')

            else:
                f.write(str(h) + '\n')

        f.close()
    
    def __compute_difference(self, i0, v0, i1, v1):

        if len(i0) == 0:
            if len(i1) == 0:
                return []
            
            else:
                return [(ChainHistory.DELTA_ADD, i1[0], v1[0])] + self.__compute_difference(i0, v0, i1[1:], v1[1:])

        elif len(i1) == 0:
            return [(ChainHistory.DELTA_REMOVE, i0[0], v0[0])] + self.__compute_difference(i0[1:], v0[1:], i1, v1)

        else:
            if (i0[0] == i1[0]):
                
                if (v0[0] == v1[0]):
                    return self.__compute_difference(i0[1:], v0[1:], i1[1:], v1[1:])
                
                else:
                    return [(ChainHistory.DELTA_CHANGE, i1[0], v1[0])] + self.__compute_difference(i0[1:], v0[1:], i1[1:], v1[1:])

            elif (i0[0] < i1[0]):

                return [(ChainHistory.DELTA_REMOVE, i0[0], v0[0])] + self.__compute_difference(i0[1:], v0[1:], i1, v1)

            else:
                return [(ChainHistory.DELTA_ADD, i1[0], v1[0])] + self.__compute_difference(i0, v0, i1[1:], v1[1:])

    def __add(self, indices, values, idx, val):

        new_indices = []
        new_values = []

        added = False
        for (i, v) in zip(indices, values):
            
            if idx == i:
                print indices, values
                print idx, val
                raise 'Index exists in add'
            elif not added and idx < i:
                new_indices.append(idx)
                new_values.append(val)
                added = True


            new_indices.append(i)
            new_values.append(v)

        if not added:
            new_indices.append(idx)
            new_values.append(val)


        return new_indices, new_values

    def __remove(self, indices, values, idx, val):

        new_indices = []
        new_values = []

        removed = False
        for (i, v) in zip(indices, values):

            if idx == i:
                removed = True
            else:
                new_indices.append(i)
                new_values.append(v)

        if not removed:
            raise 'Index missing in remove'

        return new_indices, new_values

    def __change(self, indices, values, idx, val):
    
        new_indices = []
        new_values = []

        changed = False
        for (i, v) in zip(indices, values):

            if idx == i:
                new_indices.append(i)
                new_values.append(val)
                changed = True
            else:
                new_indices.append(i)
                new_values.append(v)

        if not changed:
            raise 'Index missing in change'

        return new_indices, new_values
            
