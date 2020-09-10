#!/usr/bin/env python

#################################################################
#
#    pedigree_parser
#        Michele Berselli
#        Harvard Medical School
#        berselli.michele@gmail.com
#
#################################################################


#################################################################
#
#    LIBRARIES
#
#################################################################
import sys, os


#################################################################
#
#    Pedigree
#       -> Member
#
#################################################################
class Pedigree(object):
    ''' object to represent a pedigree '''

    class Member(object):
        ''' object to represent a member of the pedigree '''

        def __init__(self, name, sample, gender):
            ''' initialize Member object '''
            self.name = name
            self.sample = sample
            self.gender = gender
            self.parents = []
            self.children = []
        #end def __init__

        def add_parent(self, parent):
            ''' add a Member object that is parent '''
            self.parents.append(parent)
        #end def add_parent

        def add_child(self, child):
            ''' add a Member object that is child '''
            self.children.append(child)
        #end def add_child

        def is_sample(self):
            ''' '''
            if self.sample: return True
            #end if
            return False
        #end def is_sample

        def get_siblings(self):
            ''' return a list of Member objects that are siblings '''
            siblings = []
            for parent in self.parents:
                for child in parent.children:
                    if child not in siblings and child != self:
                        siblings.append(child)
                    #end if
                #end for
            #end for
            return siblings
        #end def get_siblings

        def get_parents(self):
            ''' return a list of Member objects that are parents '''
            return self.parents
        #end def get_parents

        def get_children(self):
            ''' return a list of Member objects that are children '''
            return self.children
        #end def get_children

        def get_spouses(self):
            ''' return a list of Member objects that are spouses,
            sorted by descending number of common children '''
            tmp_spouses = [] # [(num_children, spouse_obj), ...]
            for child in self.children:
                for parent in child.parents:
                    if parent != self:
                        num_children = len(self.common_children(parent))
                        if (num_children, parent) not in tmp_spouses:
                            tmp_spouses.append((num_children, parent))
                        #end if
                    #end if
                #end for
            #end for
            return [tmp_spouse[1] for tmp_spouse in sorted(tmp_spouses, key=lambda x: x[0], reverse=True)]
        #end def get_spouses

        def common_children(self, spouse):
            ''' return a list of Member objects that are children in common with spouse '''
            return list(set(self.children).intersection(set(spouse.children)))
        #end def common_children

        def has_parents(self):
            ''' check if self has two parents with samples '''
            for i, parent in enumerate(self.parents):
                if not parent.sample: # missing parent sample information
                    return False
                #end if
            #end for
            try:
                if i != 1: return False # missing one parent
                #end if
            except Exception: return False # missing both parents
            #end try
            return True
        #end def has_parents

    #end class Member

    def __init__(self, pedigree):
        ''' initialize Pedigree object,
        pedigree information must be provided as json '''
        self.members = {} # dictionary of Member objects by name
        self.samples = {} # dictionary to map sample to name
        self.parse_pedigree(pedigree)
    #end def __init__

    def add_member(self, member):
        ''' create Member object for member,
        member information must be provided as dict '''
        # Get values
        try: name = member['individual'] # individual unique ID, every member must have one
        except Exception:
            raise ValueError('\nERROR in pedigree structure, missing name for pedigree member\n')
        #end try
        try: sample = member['sample_name'] # sample ID, if available
        except Exception: sample = None
        #end try
        gender = member['gender']
        # Create object
        self.members.setdefault(name, self.Member(name, sample, gender))
        if sample: self.samples.setdefault(sample, name)
        #end if
    #end def add_member

    def add_parent(self, name, parent_name):
        ''' add parent-child relation to parent and child by name '''
        if name not in self.members:
            raise ValueError('\nERROR in pedigree structure, missing name {0} in pedigree\n'
                             .format(name))
        #end if
        if parent_name not in self.members:
            return
        #end if
        self.members[name].add_parent(self.members[parent_name])
        self.members[parent_name].add_child(self.members[name])
    #end def add_parent

    def get_member_by_sample(self, sample):
        ''' return Member object by sample '''
        try:
            return self.members[self.samples[sample]]
        except Exception:
            raise ValueError('\nERROR in pedigree structure, missing sample {0} in pedigree\n'
                             .format(sample))
        #end try
    #end def get_member_by_sample

    def get_family(self, sample):
        ''' given sample in pedigree,
        build family structure around it '''
        family = {
            'children': [],
            'parents': []
        }
        sample_obj = self.get_member_by_sample(sample)
        if sample_obj.get_children(): # sample is the center of the family,
                                      # check for parents and children
            if sample_obj.has_parents(): # parents information complete
                spouse = None
                for spouse_ in sample_obj.get_spouses():
                    if spouse_.has_parents():
                        spouse = spouse_
                        break
                    #end if
                #end for
                if not spouse:
                    raise ValueError('\nERROR in building family from pedigree, missing parents information for sample {0} spouse\n'
                            .format(sample))
                #end if
                # Create family
                family['children'] = sample_obj.common_children(spouse)
                family['parents'] = [obj for obj in sorted([sample_obj, spouse], key=lambda x: x.name)]
            else: # missing parents
                raise ValueError('\nERROR in building family from pedigree, missing parents information for sample {0}\n'
                        .format(sample))
            #end if
        else: # sample is the newest generation,
              # check for siblings, parents and grandparents
            if sample_obj.has_parents(): # parents information complete, check grandparents
                parents = sample_obj.get_parents()
                for parent in parents:
                    if not parent.has_parents():
                        raise ValueError('\nERROR in building family from pedigree, missing family information for sample {0}\n'
                                .format(sample))
                    #end if
                #end for
                # Create family
                family['children'] = parents[0].common_children(parents[1])
                family['parents'] = [obj for obj in sorted(parents, key=lambda x: x.name)]
            else:
                raise ValueError('\nERROR in building family from pedigree, missing parents information for sample {0}\n'
                        .format(sample))
            #end if
        #end if
        return family
    #end def get_family

    def parse_pedigree(self, pedigree):
        ''' read pedigree information to build Pedigree object,
        pedigree information must be provided as json '''
        # Creating Member objects
        for member in pedigree:
            self.add_member(member)
        #end for
        # Adding relations
        for member in pedigree:
            name = member['individual']
            if len(member['parents']) > 2:
                raise ValueError('\nERROR in pedigree structure, {0} has more than two parents\n'
                                 .format(name))
            #end if
            for parent_name in member['parents']:
                self.add_parent(name, parent_name)
            #end for
        #end for
    #end def parse_pedigree

#end def Pedigree
