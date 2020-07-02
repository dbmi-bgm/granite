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
            self.childrens = []
        #end def __init__

        def add_parent(self, parent):
            ''' add a Member object that is parent '''
            self.parents.append(parent)
        #end def add_parent

        def add_children(self, children):
            ''' add a Member object that is children '''
            self.childrens.append(children)
        #end def add_children

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
                for children in parent.childrens:
                    if children not in siblings and children != self:
                        siblings.append(children)
                    #end if
                #end for
            #end for
            return siblings
        #end def get_siblings

        def get_parents(self):
            ''' return a list of Member objects that are parents '''
            return self.parents
        #end def get_parents

        def get_childrens(self):
            ''' return a list of Member objects that are childrens '''
            return self.childrens
        #end def get_childrens

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
        ''' add parent-children relation to parent and children by name '''
        if name not in self.members:
            raise ValueError('\nERROR in pedigree structure, missing name {0} in pedigree\n'
                             .format(name))
        #end if
        if parent_name not in self.members:
            raise ValueError('\nERROR in pedigree structure, missing parent name {0} in pedigree\n'
                             .format(parent_name))
        #end if
        self.members[name].add_parent(self.members[parent_name])
        self.members[parent_name].add_children(self.members[name])
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
            for parent_name in member['parents']:
                self.add_parent(name, parent_name)
            #end for
        #end for
    #end def parse_pedigree

#end def Pedigree
