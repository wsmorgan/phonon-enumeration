"""Methods for testing the subroutines in the phonons module."""
import unittest as ut

def _read_output(test):
    values = []
    with open("tests/phonons/"+test) as f:
        for line in f:
            values.append(eval(line))
    return values

class TestGetArrowConcs(ut.TestCase):
    """ Tests of the get_arrow_concs subroutine."""

    def test_1(self):
        from phenum.phonons import get_arrow_concs
        params = {
            "bulk": True,
            "sizes": [],
            "lat_vecs": [],
            "nspecies": 4,
            "basis_vecs": [],
            "is_crestricted": False,
            "arrows": False,
            "concs": []
        }

        self.assertEqual(get_arrow_concs(params),[0,0,0,0])

    def test_2(self):
        from phenum.phonons import get_arrow_concs
        params = {
            "bulk": True,
            "sizes": [],
            "lat_vecs": [],
            "nspecies": 2,
            "basis_vecs": [],
            "is_crestricted": True,
            "arrows": False,
            "concs": [[1, 4, 4, 0],[2, 4, 4, 0]]
        }

        self.assertEqual(get_arrow_concs(params),[0,0])

    def test_3(self):
        from phenum.phonons import get_arrow_concs
        params = {
            "bulk": True,
            "sizes": [],
            "lat_vecs": [],
            "nspecies": 3,
            "basis_vecs": [],
            "is_crestricted": True,
            "arrows": True,
            "concs": [[0, 3, 6, 1],[3, 6, 6, 2],[0, 6, 6, 0]]
        }

        self.assertEqual(get_arrow_concs(params),[1,2,0])

    def test_4(self):
        from phenum.phonons import get_arrow_concs
        params = {
            "bulk": True,
            "sizes": [],
            "lat_vecs": [],
            "nspecies": 10,
            "basis_vecs": [],
            "is_crestricted": False,
            "arrows": False,
            "concs": []
        }

        self.assertEqual(get_arrow_concs(params),[0,0,0,0,0,0,0,0,0,0])

    def test_5(self):
        from phenum.phonons import get_arrow_concs
        params = {
            "bulk": True,
            "sizes": [],
            "lat_vecs": [],
            "nspecies": 1,
            "basis_vecs": [],
            "is_crestricted": True,
            "arrows": False,
            "concs": [[0, 4, 4, 10]]
        }

        self.assertEqual(get_arrow_concs(params),[0])

    def test_6(self):
        from phenum.phonons import get_arrow_concs
        params = {
            "bulk": True,
            "sizes": [],
            "lat_vecs": [],
            "nspecies": 5,
            "basis_vecs": [],
            "is_crestricted": True,
            "arrows": True,
            "concs": [[0, 3, 6, 1],[3, 6, 6, 2],[0, 6, 6, 0],[0, 3, 6, 10],[3, 6, 6, 1]]
        }

        self.assertEqual(get_arrow_concs(params),[1,2,0,10,1])

    def test_7(self):
        from phenum.phonons import get_arrow_concs
        params = {
            "bulk": True,
            "sizes": [],
            "lat_vecs": [],
            "nspecies": 3,
            "basis_vecs": [],
            "is_crestricted": True,
            "arrows": True,
            "concs": [[0, 3, 6, 0],[3, 6, 6, 0],[0, 6, 6, 0]]
        }

        self.assertEqual(get_arrow_concs(params),[0,0,0])

    def test_8(self):
        from phenum.phonons import get_arrow_concs
        params = {
            "bulk": True,
            "sizes": [],
            "lat_vecs": [],
            "nspecies": 4,
            "basis_vecs": [],
            "is_crestricted": False,
            "arrows": False,
            "concs": [[0, 3, 6, 3],[3, 6, 6, 2],[0, 6, 6, 1],[0, 6, 6, 1]]
        }

        self.assertEqual(get_arrow_concs(params),[0,0,0,0])

    def test_9(self):
        from phenum.phonons import get_arrow_concs
        params = {
            "bulk": True,
            "sizes": [],
            "lat_vecs": [],
            "nspecies": 1,
            "basis_vecs": [],
            "is_crestricted": True,
            "arrows": True,
            "concs": [[0, 6, 6, 5]]
        }

        self.assertEqual(get_arrow_concs(params),[5])

    def test_10(self):
        from phenum.phonons import get_arrow_concs
        params = {
            "bulk": True,
            "sizes": [],
            "lat_vecs": [],
            "nspecies": 2,
            "basis_vecs": [],
            "is_crestricted": True,
            "arrows": True,
            "concs": [[0, 3, 6, 2],[3, 6, 6, 3]]
        }

        self.assertEqual(get_arrow_concs(params),[2,3])
        
class TestArrowConcs(ut.TestCase):
    """Tests of the arrow_concs subroutine."""

    def test_1(self):
        from phenum.phonons import arrow_concs

        cList = [1, 2, 1]
        aconcs = [0, 0.4245868629437351, 0]

        self.assertEqual(arrow_concs(cList,aconcs),[[-1, 1], [-1, 3], [-1, 2], [-1, 2]])
        
    def test_2(self):
        from phenum.phonons import arrow_concs

        cList = [3]
        aconcs = [0.8205195542173467]
        out = [[-1, 1], [1, 1], [1, 1]]
        
        self.assertEqual(arrow_concs(cList,aconcs),out)
        
    def test_3(self):
        from phenum.phonons import arrow_concs

        cList = [10, 3, 1]
        aconcs = [0, 0.4989661535030203, 0]
        out = [[-1, 3], [-1, 2], [-1, 2], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1],
               [-1, 1], [-1, 1], [-1, 1], [-1, 1], [1, 2]]
        
        self.assertEqual(arrow_concs(cList,aconcs),out)
        
    def test_4(self):
        from phenum.phonons import arrow_concs

        cList = [2, 4, 1, 5, 2, 1, 1]
        aconcs = [0.9068065455664464, 0.2858477549741846, 0, 0, 0.6957735268097871, 0, 0]
        out = [[-1, 1], [-1, 3], [-1, 5], [-1, 6], [-1, 7], [-1, 2], [-1, 2], [-1, 2], [-1, 4],
               [-1, 4], [-1, 4], [-1, 4], [-1, 4], [1, 1], [1, 2], [1, 5]]
        
        self.assertEqual(arrow_concs(cList,aconcs),out)
        
    def test_5(self):
        from phenum.phonons import arrow_concs

        cList = [3]
        aconcs = [0.2871674398220775]
        out = [[-1, 1], [-1, 1], [-1, 1]]
        
        self.assertEqual(arrow_concs(cList,aconcs),out)
        
    def test_6(self):
        from phenum.phonons import arrow_concs

        cList = [3, 1]
        aconcs = [0.32514696876724436, 0]
        out = [[-1, 2], [-1, 1], [-1, 1], [-1, 1]]
        
        self.assertEqual(arrow_concs(cList,aconcs),out)
        
    def test_7(self):
        from phenum.phonons import arrow_concs

        cList = [1]
        aconcs = [0]
        out = [[-1, 1]]
        
        self.assertEqual(arrow_concs(cList,aconcs),out)
        
    def test_8(self):
        from phenum.phonons import arrow_concs

        cList = [2, 8, 3, 1]
        aconcs = [0.8244881520042212, 0.33517966472359717, 0.677253228566329, 0]
        out = [[-1, 1], [-1, 3], [-1, 4], [-1, 2], [-1, 2], [-1, 2], [-1, 2], [-1, 2], [-1, 2],
               [1, 1], [1, 2], [1, 2], [1, 3], [1, 3]]
        
        self.assertEqual(arrow_concs(cList,aconcs),out)
        
    def test_9(self):
        from phenum.phonons import arrow_concs

        cList = [4, 1, 1, 1]
        aconcs = [0, 0, 0, 0]
        out = [[-1, 2], [-1, 3], [-1, 4], [-1, 1], [-1, 1], [-1, 1], [-1, 1]]
        
        self.assertEqual(arrow_concs(cList,aconcs),out)
        
    def test_10(self):
        from phenum.phonons import arrow_concs

        cList = [18, 1]
        aconcs = [0,0]
        out = [[-1, 2], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1],
               [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1],
               [-1, 1]]
        
        self.assertEqual(arrow_concs(cList,aconcs),out)
        
    def test_11(self):
        from phenum.phonons import arrow_concs

        cList = [2, 0, 1]
        aconcs = [0,0,0]
        out = [[-1, 3], [-1, 1], [-1,1]]
        
        self.assertEqual(arrow_concs(cList,aconcs),out)
        
    def test_12(self):
        from phenum.phonons import arrow_concs

        cList = [3, 0, 2]
        aconcs = [0,0,0.5]
        out = [[-1, 3], [-1, 1], [-1, 1], [-1, 1], [1, 3]]
        
        self.assertEqual(arrow_concs(cList,aconcs),out)
        

class TestHowManyArrows(ut.TestCase):
    """Tests of the how_many_arrows subroutine."""

    def test_1(self):
        from phenum.phonons import how_many_arrows
        tcol = [[-1, 2], [-1, 2], [-1, 1], [-1, 3]]
        out = (0,0,[2,1,1])
        
        self.assertEqual(how_many_arrows(tcol),out)

    def test_2(self):
        from phenum.phonons import how_many_arrows

        tcol = [[-1, 1], [1, 1], [1, 1]]
        out = (2,1,[1,2])
        
        self.assertEqual(how_many_arrows(tcol),out)
        
    def test_3(self):
        from phenum.phonons import how_many_arrows

        tcol = [[-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1],
               [-1, 1], [-1, 2], [-1, 2], [-1, 3], [1, 2]]
        out = (1,1,[10,2,1,1])
        self.assertEqual(how_many_arrows(tcol),out)
        
    def test_4(self):
        from phenum.phonons import how_many_arrows

        tcol = [[-1, 4], [-1, 4], [-1, 4], [-1, 4], [-1, 4], [-1, 2], [-1, 2], [-1, 2], [-1, 1],
               [-1, 3], [-1, 5], [-1, 6], [-1, 7], [1, 1], [1, 2], [1, 5]]
        out = (3,3,[5,3,1,1,1,1,1,1,1,1])
        self.assertEqual(how_many_arrows(tcol),out)
        
    def test_5(self):
        from phenum.phonons import how_many_arrows

        tcol = [[-1, 1], [-1, 1], [-1, 1]]
        out = (0,0,[3])
        self.assertEqual(how_many_arrows(tcol),out)
        
    def test_6(self):
        from phenum.phonons import how_many_arrows

        tcol = [[-1, 1], [-1, 1], [-1, 1], [-1, 2]]
        out = (0,0,[3,1])
        self.assertEqual(how_many_arrows(tcol),out)
        
    def test_7(self):
        from phenum.phonons import how_many_arrows

        tcol = [[-1, 1]]
        out = (0,0,[1])
        self.assertEqual(how_many_arrows(tcol),out)
        
    def test_8(self):
        from phenum.phonons import how_many_arrows

        tcol = [[-1, 2], [-1, 2], [-1, 2], [-1, 2], [-1, 2], [-1, 2], [-1, 1], [-1, 3], [-1, 4],
               [1, 2], [1, 2], [1, 3], [1, 3], [1, 1]]
        out = (5,3,[6,1,1,1,2,2,1])
        self.assertEqual(how_many_arrows(tcol),out)
        
    def test_9(self):
        from phenum.phonons import how_many_arrows

        tcol = [[-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 2], [-1, 3], [-1, 4]]
        out = (0,0,[4,1,1,1])
        self.assertEqual(how_many_arrows(tcol),out)
        
    def test_10(self):
        from phenum.phonons import how_many_arrows

        tcol = [[-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1],
                [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1],
                [-1, 2]]
        out = (0,0,[18,1])
        self.assertEqual(how_many_arrows(tcol),out)

class TestEnumSys(ut.TestCase):
    """Tests of the enum_sys subroutine."""
    
    def test_1(self):
        from phenum.phonons import enum_sys
        from numpy import array
        groupfile = "tests/phonons/test_group.1"
        concs = [1,2]
        a_cons = [0,0]
        num_wanted = 1
        HNF = array([1,0,1,0,2,3])
        params ={'bulk': True, 'nspecies': 2, 'concs': [], 'basis_vecs': [[0.0, 0.0, 0.0]], 'sizes': [1, 11], 'lat_vecs': [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]], 'arrows': False, 'is_crestricted': False}
        out = [[[-1, 1], [-1, 2], [-1, 2]]]
        self.assertEqual(enum_sys(groupfile,concs,a_cons,num_wanted,HNF,params,True),out)

    def test_2(self):
        from phenum.phonons import enum_sys
        from numpy import array
        groupfile = "tests/phonons/test_group.2"
        concs = [3,3]
        a_cons = [0,0]
        num_wanted = 3
        HNF = array([1,0,1,0,0,6])
        params = {'bulk': True, 'nspecies': 2, 'concs': [], 'basis_vecs': [[0.0, 0.0, 0.0]], 'sizes': [1, 11], 'lat_vecs': [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]], 'arrows': False, 'is_crestricted': False}
        out = [[[-1, 1], [-1, 1], [-1, 1], [-1, 2], [-1, 2], [-1, 2]], [[-1, 1], [-1, 1], [-1, 2], [-1, 1], [-1, 2], [-1, 2]], [[-1, 1], [-1, 2], [-1, 1], [-1, 2], [-1, 1], [-1, 2]]]
        self.assertEqual(enum_sys(groupfile,concs,a_cons,num_wanted,HNF,params,True),out)

    def test_3(self):
        from numpy import array
        from phenum.phonons import enum_sys
        groupfile = "tests/phonons/test_group.3"
        concs = [4,3]
        a_cons = [0,0]
        num_wanted = 4
        HNF = array([1,0,1,1,2,7])
        params = {'bulk': True, 'nspecies': 2, 'concs': [], 'basis_vecs': [[0.0, 0.0, 0.0]], 'sizes': [1, 11], 'lat_vecs': [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]], 'arrows': False, 'is_crestricted': False}
        out = _read_output("enum_sys.out.3")
                
        self.assertEqual(enum_sys(groupfile,concs,a_cons,num_wanted,HNF,params,True),out)

    def test_4(self):
        from phenum.phonons import enum_sys
        from numpy import array
        groupfile = "tests/phonons/test_group.4"
        concs = [3,4]
        a_cons = [0,0]
        num_wanted = 2
        HNF = array([1,0,1,1,3,7])
        params = {'bulk': True, 'nspecies': 2, 'concs': [], 'basis_vecs': [[0.0, 0.0, 0.0]], 'sizes': [1, 11], 'lat_vecs': [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]], 'arrows': False, 'is_crestricted': False}
        out = [[[-1, 1], [-1, 1], [-1, 1], [-1, 2], [-1, 2], [-1, 2], [-1, 2]], [[-1, 1], [-1, 1], [-1, 2], [-1, 1], [-1, 2], [-1, 2], [-1, 2]]]
        self.assertEqual(enum_sys(groupfile,concs,a_cons,num_wanted,HNF,params,True),out)

    def test_5(self):
        from phenum.phonons import enum_sys
        from numpy import array
        groupfile = "tests/phonons/test_group.5"
        concs = [4,4]
        a_cons = [0,0]
        num_wanted = 10
        HNF = array([1,0,2,0,0,4])
        params = {'bulk': True, 'nspecies': 2, 'concs': [], 'basis_vecs': [[0.0, 0.0, 0.0]], 'sizes': [1, 11], 'lat_vecs': [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]], 'arrows': False, 'is_crestricted': False}
        out = [[[-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 2], [-1, 2], [-1, 2], [-1, 2]], [[-1, 1], [-1, 1], [-1, 1], [-1, 2], [-1, 1], [-1, 2], [-1, 2], [-1, 2]], [[-1, 1], [-1, 1], [-1, 1], [-1, 2], [-1, 2], [-1, 1], [-1, 2], [-1, 2]], [[-1, 1], [-1, 1], [-1, 1], [-1, 2], [-1, 2], [-1, 2], [-1, 2], [-1, 1]], [[-1, 1], [-1, 1], [-1, 2], [-1, 2], [-1, 1], [-1, 1], [-1, 2], [-1, 2]], [[-1, 1], [-1, 1], [-1, 2], [-1, 2], [-1, 1], [-1, 2], [-1, 1], [-1, 2]], [[-1, 1], [-1, 1], [-1, 2], [-1, 2], [-1, 1], [-1, 2], [-1, 2], [-1, 1]], [[-1, 1], [-1, 1], [-1, 2], [-1, 2], [-1, 2], [-1, 2], [-1, 1], [-1, 1]], [[-1, 1], [-1, 2], [-1, 1], [-1, 2], [-1, 1], [-1, 2], [-1, 1], [-1, 2]], [[-1, 1], [-1, 2], [-1, 1], [-1, 2], [-1, 2], [-1, 1], [-1, 2], [-1, 1]]]
        self.assertEqual(enum_sys(groupfile,concs,a_cons,num_wanted,HNF,params,True),out)

    def test_6(self):
        from phenum.phonons import enum_sys
        from numpy import array
        groupfile = None
        concs = [3,1,2]
        a_cons = [0.0,0.5,0.25]
        num_wanted = 6
        HNF = array([1,0,1,0,0,6])
        params = {'bulk': True, 'nspecies': 3, 'concs': [[1.0, 6.0, 12.0, 0.0], [1.0, 9.0, 12.0, 0.5], [1.0, 12.0, 12.0, 0.25]], 'basis_vecs': [[0.0, 0.0, 0.0]], 'sizes': [6, 6], 'lat_vecs': [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]], 'arrows': True, 'is_crestricted': True}
        out = _read_output("enum_sys.out.6")

        self.assertEqual(enum_sys(groupfile,concs,a_cons,num_wanted,HNF,params,True),out)

    # def test_7(self):
    #     from phenum.phonons import enum_sys
    #     from numpy import array
    #     groupfile = None
    #     concs = [1,4,1]
    #     a_cons = [0.0,0.5,0.25]
    #     num_wanted = 124
    #     HNF = array([1,0,1,0,5,6])
    #     params = {'bulk': True, 'nspecies': 3, 'concs': [[1.0, 6.0, 12.0, 0.0], [1.0, 9.0, 12.0, 0.5], [1.0, 12.0, 12.0, 0.25]], 'basis_vecs': [[0.0, 0.0, 0.0]], 'sizes': [3, 3], 'lat_vecs': [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]], 'arrows': True, 'is_crestricted': True}
    #     out = _read_output("enum_sys.out.7")

    #     self.assertEqual(enum_sys(groupfile,concs,a_cons,num_wanted,HNF,params,True),out)

    # def test_8(self):
    #     from phenum.phonons import enum_sys
    #     from numpy import array
    #     groupfile = None
    #     concs = [1,3]
    #     a_cons = [0.0,0.75]
    #     num_wanted = 19
    #     HNF = array([1,0,1,0,1,2])
    #     params = {'bulk': True, 'nspecies': 2, 'concs': [[1.0, 6.0, 12.0, 0.0], [1.0, 9.0, 12.0, 0.75]], 'basis_vecs': [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], 'sizes': [2, 2], 'lat_vecs': [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], 'arrows': True, 'is_crestricted': True}
    #     out = _read_output("enum_sys.out.8")
    #     self.assertEqual(enum_sys(groupfile,concs,a_cons,num_wanted,HNF,params,True),out)

    # def test_9(self):
    #     from phenum.phonons import enum_sys
    #     from numpy import array
    #     groupfile = None
    #     concs = [2,5]
    #     a_cons = [0.25,0.75]
    #     num_wanted = 738
    #     HNF = array([1,0,1,0,0,7])
    #     params = {'bulk': True, 'nspecies': 2, 'concs': [[1.0, 6.0, 12.0, 0.25], [1.0, 9.0, 12.0, 0.75]], 'basis_vecs': [[0.0, 0.0, 0.0]], 'sizes': [7, 7], 'lat_vecs': [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], 'arrows': True, 'is_crestricted': True}
    #     out = _read_output("enum_sys.out.9")
    #     self.assertEqual(enum_sys(groupfile,concs,a_cons,num_wanted,HNF,params,True),out)

    # def test_10(self):
    #     from phenum.phonons import enum_sys
    #     from numpy import array
    #     groupfile = None
    #     concs = [1,1,1,1]
    #     a_cons = [0.0,0.0,1.0,1.0]
    #     num_wanted = 36
    #     HNF = array([1,0,2,0,0,2])
    #     params = {'bulk': True, 'nspecies': 4, 'concs': [[0.0, 4.0, 4.0, 0.0], [0.0, 4.0, 4.0, 0.0], [0.0, 4.0, 4.0, 1.0], [0.0, 4.0, 4.0, 1.0]], 'basis_vecs': [[0.0, 0.0, 0.0]], 'sizes': [4, 4], 'lat_vecs': [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], 'arrows': True, 'is_crestricted': True}
    #     out = _read_output("enum_sys.out.10")
    #     self.assertEqual(enum_sys(groupfile,concs,a_cons,num_wanted,HNF,params,True),out)

    # def test_11(self):
    #     from phenum.phonons import enum_sys
    #     from numpy import array
    #     groupfile = None
    #     concs = [1,0,1,1,1]
    #     a_cons = [0.0,0.0,0.0,1.0,1.0]
    #     num_wanted = 36
    #     HNF = array([1,0,2,0,0,2])
    #     params = {'bulk': True, 'nspecies': 4, 'concs': [[0.0, 4.0, 4.0, 0.0], [0.0, 4.0, 4.0, 0.0], [0.0, 4.0, 4.0, 1.0], [0.0, 4.0, 4.0, 1.0]], 'basis_vecs': [[0.0, 0.0, 0.0]], 'sizes': [4, 4], 'lat_vecs': [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], 'arrows': True, 'is_crestricted': True}
    #     out = _read_output("enum_sys.out.11")

    #     self.assertEqual(enum_sys(groupfile,concs,a_cons,num_wanted,HNF,params,True),out)

    # def test_12(self):
    #     from phenum.phonons import enum_sys
    #     from numpy import array
    #     groupfile = None
    #     concs = [2,0,5]
    #     a_cons = [0.25,0.0,0.75]
    #     num_wanted = 738
    #     HNF = array([1,0,1,0,0,7])
    #     params = {'bulk': True, 'nspecies': 2, 'concs': [[1.0, 6.0, 12.0, 0.25], [1.0, 9.0, 12.0, 0.75]], 'basis_vecs': [[0.0, 0.0, 0.0]], 'sizes': [7, 7], 'lat_vecs': [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], 'arrows': True, 'is_crestricted': True}
    #     out = _read_output("enum_sys.out.12")

    #     self.assertEqual(enum_sys(groupfile,concs,a_cons,num_wanted,HNF,params,True),out)

    # def test_13(self):
    #     from phenum.phonons import enum_sys
    #     from numpy import array
    #     groupfile = None
    #     concs = [0,0,7]
    #     a_cons = [0.0,0.0,0.0]
    #     num_wanted = 738
    #     HNF = array([1,0,1,0,0,7])
    #     params = {'bulk': True, 'nspecies': 2, 'concs': [[1.0, 6.0, 12.0, 0.25], [1.0, 9.0, 12.0, 0.75]], 'basis_vecs': [[0.0, 0.0, 0.0]], 'sizes': [7, 7], 'lat_vecs': [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], 'arrows': True, 'is_crestricted': True}
    #     out = []

    #     self.assertEqual(enum_sys(groupfile,concs,a_cons,num_wanted,HNF,params,False),out)

    # def test_14(self):
    #     from phenum.phonons import enum_sys
    #     from numpy import array
    #     groupfile = None
    #     concs = [1,0,1,1,1]
    #     a_cons = [0.0,0.0,0.0,1.0,1.0]
    #     num_wanted = 36
    #     HNF = array([1,0,2,0,0,2])
    #     params = {'bulk': True, 'nspecies': 4, 'concs': [[0.0, 4.0, 4.0, 0.0], [0.0, 4.0, 4.0, 0.0], [0.0, 4.0, 4.0, 1.0], [0.0, 4.0, 4.0, 1.0]], 'basis_vecs': [[0.0, 0.0, 0.0]], 'sizes': [4, 4], 'lat_vecs': [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], 'arrows': True, 'is_crestricted': True}
    #     out = _read_output("enum_sys.out.14")

    #     self.assertEqual(enum_sys(groupfile,concs,a_cons,num_wanted,HNF,params,False),out)

    def test_15(self):
        from phenum.phonons import enum_sys
        from numpy import array
        groupfile = "tests/phonons/test_group.5"
        concs = [4,4]
        a_cons = [0,0]
        num_wanted = 10
        HNF = array([1,0,2,0,0,4])
        params = {'bulk': True, 'nspecies': 2, 'concs': [], 'basis_vecs': [[0.0, 0.0, 0.0]], 'sizes': [1, 11], 'lat_vecs': [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]], 'arrows': False, 'is_crestricted': False}
        out = _read_output("enum_sys.out.15")

        self.assertEqual(enum_sys(groupfile,concs,a_cons,num_wanted,HNF,params,False),out)

class TestAddArrows(ut.TestCase):
    """Tests of the add_arrows subroutine."""

    def test_1(self):
        from phenum.grouptheory import get_sym_group
        from phenum.phonons import add_arrows
        col = [[-1, 2], [-1, 1], [1, 2], [-1, 2], [1, 2], [-1, 3]]
        agroup = _read_output("add_arrow_group.in.1")
        dim = 6
        out = [[[-1, 2], [-1, 1], [0, 2], [-1, 2], [0, 2], [-1, 3]], [[-1, 2], [-1, 1], [1, 2], [-1, 2], [0, 2], [-1, 3]], [[-1, 2], [-1, 1], [2, 2], [-1, 2], [0, 2], [-1, 3]], [[-1, 2], [-1, 1], [5, 2], [-1, 2], [0, 2], [-1, 3]], [[-1, 2], [-1, 1], [0, 2], [-1, 2], [2, 2], [-1, 3]], [[-1, 2], [-1, 1], [2, 2], [-1, 2], [2, 2], [-1, 3]], [[-1, 2], [-1, 1], [3, 2], [-1, 2], [2, 2], [-1, 3]], [[-1, 2], [-1, 1], [4, 2], [-1, 2], [2, 2], [-1, 3]]]
        self.assertEqual(add_arrows(col,agroup,dim,agroup[0:len(col)],supers=True),out)

    def test_2(self):
        from phenum.grouptheory import get_sym_group
        from phenum.phonons import add_arrows
        col = [[-1, 1], [1, 2], [-1, 2], [1, 2]]
        agroup = _read_output("add_arrow_group.in.2")
        dim = 6
        out = [[[-1, 1], [0, 2], [-1, 2], [0, 2]], [[-1, 1], [1, 2], [-1, 2], [0, 2]], [[-1, 1], [5, 2], [-1, 2], [0, 2]], [[-1, 1], [0, 2], [-1, 2], [1, 2]], [[-1, 1], [1, 2], [-1, 2], [1, 2]], [[-1, 1], [2, 2], [-1, 2], [1, 2]], [[-1, 1], [4, 2], [-1, 2], [1, 2]], [[-1, 1], [5, 2], [-1, 2], [1, 2]], [[-1, 1], [0, 2], [-1, 2], [5, 2]], [[-1, 1], [1, 2], [-1, 2], [5, 2]], [[-1, 1], [5, 2], [-1, 2], [5, 2]]]
        self.assertEqual(add_arrows(col,agroup,dim,agroup[0:len(col)],supers=True),out)

    def test_3(self):
        from phenum.grouptheory import get_sym_group
        from phenum.phonons import add_arrows
        col = [[-1, 1], [1, 2], [1, 2], [-1, 1], [-1, 2], [1, 2], [-1, 2]]
        agroup = _read_output("add_arrow_group.in.3")
        dim = 6
        out = [[[-1, 1], [0, 2], [0, 2], [-1, 1], [-1, 2], [0, 2], [-1, 2]], [[-1, 1], [1, 2], [0, 2], [-1, 1], [-1, 2], [0, 2], [-1, 2]], [[-1, 1], [5, 2], [0, 2], [-1, 1], [-1, 2], [0, 2], [-1, 2]], [[-1, 1], [0, 2], [1, 2], [-1, 1], [-1, 2], [0, 2], [-1, 2]], [[-1, 1], [1, 2], [1, 2], [-1, 1], [-1, 2], [0, 2], [-1, 2]], [[-1, 1], [2, 2], [1, 2], [-1, 1], [-1, 2], [0, 2], [-1, 2]], [[-1, 1], [4, 2], [1, 2], [-1, 1], [-1, 2], [0, 2], [-1, 2]], [[-1, 1], [5, 2], [1, 2], [-1, 1], [-1, 2], [0, 2], [-1, 2]], [[-1, 1], [0, 2], [5, 2], [-1, 1], [-1, 2], [0, 2], [-1, 2]], [[-1, 1], [1, 2], [5, 2], [-1, 1], [-1, 2], [0, 2], [-1, 2]], [[-1, 1], [5, 2], [5, 2], [-1, 1], [-1, 2], [0, 2], [-1, 2]], [[-1, 1], [0, 2], [0, 2], [-1, 1], [-1, 2], [1, 2], [-1, 2]], [[-1, 1], [1, 2], [0, 2], [-1, 1], [-1, 2], [1, 2], [-1, 2]], [[-1, 1], [2, 2], [0, 2], [-1, 1], [-1, 2], [1, 2], [-1, 2]], [[-1, 1], [4, 2], [0, 2], [-1, 1], [-1, 2], [1, 2], [-1, 2]], [[-1, 1], [5, 2], [0, 2], [-1, 1], [-1, 2], [1, 2], [-1, 2]], [[-1, 1], [0, 2], [1, 2], [-1, 1], [-1, 2], [1, 2], [-1, 2]], [[-1, 1], [1, 2], [1, 2], [-1, 1], [-1, 2], [1, 2], [-1, 2]], [[-1, 1], [2, 2], [1, 2], [-1, 1], [-1, 2], [1, 2], [-1, 2]], [[-1, 1], [4, 2], [1, 2], [-1, 1], [-1, 2], [1, 2], [-1, 2]], [[-1, 1], [0, 2], [2, 2], [-1, 1], [-1, 2], [1, 2], [-1, 2]], [[-1, 1], [2, 2], [2, 2], [-1, 1], [-1, 2], [1, 2], [-1, 2]], [[-1, 1], [3, 2], [2, 2], [-1, 1], [-1, 2], [1, 2], [-1, 2]], [[-1, 1], [4, 2], [2, 2], [-1, 1], [-1, 2], [1, 2], [-1, 2]], [[-1, 1], [0, 2], [4, 2], [-1, 1], [-1, 2], [1, 2], [-1, 2]], [[-1, 1], [4, 2], [4, 2], [-1, 1], [-1, 2], [1, 2], [-1, 2]], [[-1, 1], [0, 2], [5, 2], [-1, 1], [-1, 2], [1, 2], [-1, 2]]]
        self.assertEqual(add_arrows(col,agroup,dim,agroup[0:len(col)],supers=True),out)

    def test_4(self):
        from phenum.grouptheory import get_sym_group
        from phenum.phonons import add_arrows
        col = [[-1, 1], [1, 3], [1, 4], [-1, 2]]
        agroup = _read_output("add_arrow_group.in.4")
        dim = 6
        out = [[[-1, 1], [0, 3], [0, 4], [-1, 2]], [[-1, 1], [1, 3], [0, 4], [-1, 2]], [[-1, 1], [2, 3], [0, 4], [-1, 2]], [[-1, 1], [5, 3], [0, 4], [-1, 2]], [[-1, 1], [0, 3], [1, 4], [-1, 2]], [[-1, 1], [1, 3], [1, 4], [-1, 2]], [[-1, 1], [2, 3], [1, 4], [-1, 2]], [[-1, 1], [4, 3], [1, 4], [-1, 2]], [[-1, 1], [0, 3], [2, 4], [-1, 2]], [[-1, 1], [1, 3], [2, 4], [-1, 2]], [[-1, 1], [2, 3], [2, 4], [-1, 2]], [[-1, 1], [3, 3], [2, 4], [-1, 2]]]
        self.assertEqual(add_arrows(col,agroup,dim,agroup[0:len(col)],supers=True),out)

    def test_5(self):
        from phenum.grouptheory import get_sym_group
        from phenum.phonons import add_arrows
        col = [[-1, 1], [-1, 2], [1, 2], [-1, 1], [-1, 1], [-1, 1]]
        agroup = _read_output("add_arrow_group.in.5")
        dim = 6
        out = [[[-1, 1], [-1, 2], [0, 2], [-1, 1], [-1, 1], [-1, 1]], [[-1, 1], [-1, 2], [1, 2], [-1, 1], [-1, 1], [-1, 1]], [[-1, 1], [-1, 2], [3, 2], [-1, 1], [-1, 1], [-1, 1]], [[-1, 1], [-1, 2], [5, 2], [-1, 1], [-1, 1], [-1, 1]]]
        self.assertEqual(add_arrows(col,agroup,dim,agroup[0:len(col)],supers=True),out)

    def test_6(self):
        from phenum.grouptheory import get_sym_group
        from phenum.phonons import add_arrows
        col = [[-1, 1], [-1, 1], [-1, 2], [1, 2], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1]]
        agroup = _read_output("add_arrow_group.in.6")
        dim = 6
        out = [[[-1, 1], [-1, 1], [-1, 2], [0, 2], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1]], [[-1, 1], [-1, 1], [-1, 2], [1, 2], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1]], [[-1, 1], [-1, 1], [-1, 2], [3, 2], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1]], [[-1, 1], [-1, 1], [-1, 2], [5, 2], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1]]]
        self.assertEqual(add_arrows(col,agroup,dim,agroup[0:len(col)],supers=True),out)

    def test_7(self):
        from phenum.grouptheory import get_sym_group
        from phenum.phonons import add_arrows
        col = [[-1, 3], [-1, 3], [-1, 1], [-1, 3], [-1, 2], [-1, 3], [-1, 3], [-1, 1], [-1, 3], [1, 2]]
        agroup = _read_output("add_arrow_group.in.7")
        dim = 6
        out = [[[-1, 3], [-1, 3], [-1, 1], [-1, 3], [-1, 2], [-1, 3], [-1, 3], [-1, 1], [-1, 3], [0, 2]], [[-1, 3], [-1, 3], [-1, 1], [-1, 3], [-1, 2], [-1, 3], [-1, 3], [-1, 1], [-1, 3], [1, 2]], [[-1, 3], [-1, 3], [-1, 1], [-1, 3], [-1, 2], [-1, 3], [-1, 3], [-1, 1], [-1, 3], [5, 2]]]
        self.assertEqual(add_arrows(col,agroup,dim,agroup[0:len(col)],supers=True),out)

    def test_8(self):
        from phenum.grouptheory import get_sym_group
        from phenum.phonons import add_arrows
        col = [[-1, 3], [-1, 3], [-1, 1], [-1, 3], [-1, 2], [-1, 3], [-1, 3], [1, 2], [-1, 3], [-1, 1]]
        agroup = _read_output("add_arrow_group.in.8")
        dim = 6
        out = [[[-1, 3], [-1, 3], [-1, 1], [-1, 3], [-1, 2], [-1, 3], [-1, 3], [0, 2], [-1, 3], [-1, 1]], [[-1, 3], [-1, 3], [-1, 1], [-1, 3], [-1, 2], [-1, 3], [-1, 3], [1, 2], [-1, 3], [-1, 1]], [[-1, 3], [-1, 3], [-1, 1], [-1, 3], [-1, 2], [-1, 3], [-1, 3], [5, 2], [-1, 3], [-1, 1]]]
        self.assertEqual(add_arrows(col,agroup,dim,agroup[0:len(col)],supers=True),out)

    def test_9(self):
        from phenum.grouptheory import get_sym_group
        from phenum.phonons import add_arrows
        col = [[-1, 1], [-1, 3], [-1, 2], [1, 4], [1, 2], [-1, 4], [-1, 3], [-1, 1]]
        agroup = _read_output("add_arrow_group.in.9")
        dim = 6
        out = [[[-1, 1], [-1, 3], [-1, 2], [0, 4], [0, 2], [-1, 4], [-1, 3], [-1, 1]], [[-1, 1], [-1, 3], [-1, 2], [1, 4], [0, 2], [-1, 4], [-1, 3], [-1, 1]], [[-1, 1], [-1, 3], [-1, 2], [2, 4], [0, 2], [-1, 4], [-1, 3], [-1, 1]], [[-1, 1], [-1, 3], [-1, 2], [5, 4], [0, 2], [-1, 4], [-1, 3], [-1, 1]], [[-1, 1], [-1, 3], [-1, 2], [0, 4], [1, 2], [-1, 4], [-1, 3], [-1, 1]], [[-1, 1], [-1, 3], [-1, 2], [1, 4], [1, 2], [-1, 4], [-1, 3], [-1, 1]], [[-1, 1], [-1, 3], [-1, 2], [2, 4], [1, 2], [-1, 4], [-1, 3], [-1, 1]], [[-1, 1], [-1, 3], [-1, 2], [4, 4], [1, 2], [-1, 4], [-1, 3], [-1, 1]], [[-1, 1], [-1, 3], [-1, 2], [0, 4], [2, 2], [-1, 4], [-1, 3], [-1, 1]], [[-1, 1], [-1, 3], [-1, 2], [1, 4], [2, 2], [-1, 4], [-1, 3], [-1, 1]], [[-1, 1], [-1, 3], [-1, 2], [2, 4], [2, 2], [-1, 4], [-1, 3], [-1, 1]], [[-1, 1], [-1, 3], [-1, 2], [3, 4], [2, 2], [-1, 4], [-1, 3], [-1, 1]]]
        self.assertEqual(add_arrows(col,agroup,dim,agroup[0:len(col)],supers=True),out)

    def test_10(self):
        from phenum.grouptheory import get_sym_group
        from phenum.phonons import add_arrows
        col = [[-1, 1], [-1, 2], [-1, 3], [1, 2], [-1, 3], [-1, 4], [-1, 1], [1, 4]]
        agroup = _read_output("add_arrow_group.in.10")
        dim = 6
        out = [[[-1, 1], [-1, 2], [-1, 3], [0, 2], [-1, 3], [-1, 4], [-1, 1], [0, 4]], [[-1, 1], [-1, 2], [-1, 3], [1, 2], [-1, 3], [-1, 4], [-1, 1], [0, 4]], [[-1, 1], [-1, 2], [-1, 3], [2, 2], [-1, 3], [-1, 4], [-1, 1], [0, 4]], [[-1, 1], [-1, 2], [-1, 3], [5, 2], [-1, 3], [-1, 4], [-1, 1], [0, 4]], [[-1, 1], [-1, 2], [-1, 3], [0, 2], [-1, 3], [-1, 4], [-1, 1], [1, 4]], [[-1, 1], [-1, 2], [-1, 3], [1, 2], [-1, 3], [-1, 4], [-1, 1], [1, 4]], [[-1, 1], [-1, 2], [-1, 3], [2, 2], [-1, 3], [-1, 4], [-1, 1], [1, 4]], [[-1, 1], [-1, 2], [-1, 3], [4, 2], [-1, 3], [-1, 4], [-1, 1], [1, 4]], [[-1, 1], [-1, 2], [-1, 3], [5, 2], [-1, 3], [-1, 4], [-1, 1], [1, 4]], [[-1, 1], [-1, 2], [-1, 3], [0, 2], [-1, 3], [-1, 4], [-1, 1], [2, 4]], [[-1, 1], [-1, 2], [-1, 3], [1, 2], [-1, 3], [-1, 4], [-1, 1], [2, 4]], [[-1, 1], [-1, 2], [-1, 3], [2, 2], [-1, 3], [-1, 4], [-1, 1], [2, 4]], [[-1, 1], [-1, 2], [-1, 3], [3, 2], [-1, 3], [-1, 4], [-1, 1], [2, 4]], [[-1, 1], [-1, 2], [-1, 3], [5, 2], [-1, 3], [-1, 4], [-1, 1], [2, 4]], [[-1, 1], [-1, 2], [-1, 3], [0, 2], [-1, 3], [-1, 4], [-1, 1], [5, 4]], [[-1, 1], [-1, 2], [-1, 3], [1, 2], [-1, 3], [-1, 4], [-1, 1], [5, 4]], [[-1, 1], [-1, 2], [-1, 3], [2, 2], [-1, 3], [-1, 4], [-1, 1], [5, 4]], [[-1, 1], [-1, 2], [-1, 3], [5, 2], [-1, 3], [-1, 4], [-1, 1], [5, 4]]]
        self.assertEqual(add_arrows(col,agroup,dim,agroup[0:len(col)],supers=True),out)
