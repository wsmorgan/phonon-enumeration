#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
AUTHORS: Conrad W. Rosenbrock, Wiley S. Morgan (June 2015)

Classes to support the calculation of coefficients for specific terms in a product
of multinomials. Construct a product class by specifying the exponent and target term
and then add multinomials using the Product instance's append(). The coefficient is
then available from the coeff().
"""
from functools import reduce
class Sequence(object):
    """Represents an exponent-limited sequence with a single root. Here,
    sequence represents a sequence of integer values k_1, k_2...k_j
    that are the exponents of a single term in a multinomial. The root
    of the sequence is one of the k_i; its children become sets of
    sequences including variables to the right of i.

    Attributes:
        used (int): The sum of the exponents that have already been used.
        parent (Sequence): A Sequence instance for the variable previouse to this one.
        kids (Sequence): The next possible Sequence instances.
        varcount (int): The number of variables in the polynomial.
    """
    def __init__(self, root, possibles, i, powersum, targets, parent=None):
        """Initializes a sequence collector for a variable. 'Term' refers to a product
        of variables like x^i.y^j.z^k, then x has index 0, y has 1, etc.

        Args:
            root (int): The exponent of the variable to the left in the multinomial term.
            possibles (list): A list of possible values for each variable in the multinomial.
            i (int): the index of the variable being sequenced in the term.
            powersum (int): the maximum value that the sum of exponents in the sequence 
              is allowed to have.

            parent (Sequence): A Sequence instance for the variable the *left* of this one
              (i.e. has index i-1).
        """
        self._root = root
        self.used = root + (0 if parent is None else parent.used)
        self.parent = parent

        #We only keep recursively defining sequences until we run out of variables in
        #the term. Possibles is a list of possible exponents for each variable in the
        #term and has the same number of items as variables in the term.
        if i < len(targets):
            #Filter the possible values for the variable being considered based on the
            #exponent of the multinomial. When multinomials are expanded, the sum of
            #the exponents in any term must be less than the exponent on the multinomial
            #times the maximum power of any of its (unexpanded) terms.

            #We find all the possible values for this variable by ensuring that:
            # 1) it's exponent is compatible with the exponents of all variables to the left of it.
            # 2) the exponent we are suggesting is in the list of possible values for the variable.
            # 3) the exponent remains positive.
            self.kids = [Sequence(p-root, possibles, i+1, powersum, targets, self) 
                         for p in possibles if p-root >= 0
                         and p-root <= targets[i]
                         and abs(p - root) <= powersum-self.used 
                         and abs(p-self.used) % possibles[1] == 0]
        else:
            self.kids = []

        self._kidcount = None
        self.varcount = len(targets)

    @property
    def kidcount(self): #pragma: no cover
        """Returns the number of children and grandchildren to the last generation.

        Returns:
            _kidcount (int): The number of descendents of the parent.
        """
        if self._kidcount is None:
            _kidcount = sum([k.kidcount for k in self.kids])
            if _kidcount == 0:
                _kidcount = len(self.kids)
        return _kidcount

    def expand(self, depth=0):
        """Recursively generates a list of all relevant sequences for this multinomial term.

        Args:
            depth (int, optional): The current depth in the recursion.

        Returns:
            sequences (list): The relevant sequences for the multinomial.
        """
        #Iterate through the child sequences and add their variable root values if
        #the total sequence sums to the target.
        sequences = []
        for kid in self.kids:
            for seq in kid.expand(depth+1):
                #Here is where the recursion happens; we add the sequence of this variable's
                #children to the right of this root.
                sequences.append((self._root,) + seq)

        if len(self.kids) == 0:
            if depth == self.varcount-1:
                return [(self._root,)]
            else: #pragma: no cover
                return [(self._root,) + (0,)*(self.varcount-(depth+1))]
        else:
            return sequences

    def expand_noappend(self, sequences, start, varindex): #pragma: no cover
        """Implements an expansion that doesn't use python's append.

        Args:
            sequences (list): A list of the relevant sequences for the multinomials.
            start (int): Starting index.
            varindex (int): The current variabl index.

        Returns:
            sequences (list): The relevant sequences for the multinomial.

        Raises:
            ValueError: if self.kidcount is zero.
        """
        if len(sequences) == 0:
            if self.kidcount == 0:
                raise ValueError("This can't happen!")
            sequences = [[0]*self.varcount for i in range(self.kidcount)]

        cursor = start
        for kid in self.kids:
            kid.expand_noappend(sequences, cursor, varindex+1)
            cursor += kid.kidcount

            if varindex == self.varcount-1:
                cursor += 1

            for k in range(start, cursor):
                sequences[k][varindex-1] = self._root
            start = cursor

        if self.kidcount == 0:
            sequences[cursor][varindex-1] = self._root
            
        return sequences

class Product(object):
    """Represents a product of multinomials for which only a single term is interesting.

    Attributes:
        coefficient (int): The scalar integer multiplying this product of multinomials.
        targets (list): A list of exponents for the terms in the products.
        multinoms (list): A list of multinomials.
    """
    def __init__(self, coefficient, targets):
        """Initializes the empty product of multinomials.

        Args:
            coefficient (int): The scalar integer multiplying this product of multinomials.
            targets (list): A list of exponents for the only interesting term in the 
              product. The list is in the order that the variables appear in each multinomial.
        """
        self.coefficient = coefficient
        self.targets = targets
        self.multinoms = []

    def coeff(self):
        """Returns the coefficient of the term with the target exponents if all the multinomials
        in the product were expanded and had their terms collected.

        Returns:
            coefficient (int): The coefficient of th term with the target exponents.
        """
        #If this is an isolated multinomial, we only need to check the coefficient of the target
        #term.
        if len(self.multinoms) == 1:
            if all(self.multinoms[0].power-t>0 for t in self.targets):
                return 0
            else:
                return self.multinoms[0].nchoosekm(self.targets)*self.coefficient
        
        from itertools import product
        #Get a list of the possible exponents for each variable in each of the multinomials.
        #We start with the first variable and choose only those combinations of exponents
        #across *all* the multinomials that give the correct target exponent for that variable.
        possibles = [n.possible_powers for n in self.multinoms]
        
        seq0 = [s for s in product(*possibles) if sum(s) == self.targets[0]]
        #Next, we construct Sequence instances for each of the first variable compatible
        #possibilities and follow them through to the other variables.
        coeffs = 0
        for seq in seq0:
            mnseq = []
            #Each sequence calculated from the first variable has an entry for each multinomial
            #in this product. The Sequence instances construct smart sequences for the remaining
            #variables in each multinomial separately
            len_seq = len(seq)
            for i in range(len_seq):
                varseq = Sequence(seq[i], possibles[i], 1, self.multinoms[i].powersum, self.targets)
                mnseq.append(varseq.expand())
            coeffs += self._sum_sequences(mnseq)

        return int(coeffs)*self.coefficient

    def _sum_sequences(self, mnseq):
        """Sums all the possible combinations of relevant sequences based of the variable 
        sequence lists specified.
        
        Args:
            mnseq (list): a list of possible variable sequences in each multinomial (one for 
              each multinomial) that might contribute to the correct target variable.

        Returns:
            coeffs (int): The coefficients of the sequence.
        """
        from itertools import product
        from operator import mul
        #We can also filter the sequences by enforcing the constraints that the exponents
        #correctly reproduce the target across all the multionmials in the product. Get hold
        #of all the combinations of sequences across the multinomials and check each for
        #conformance to the targets.
        coeffs = 0
        for seq in product(*mnseq):
            expsum = [sum(zs) for zs in zip(*seq)]
            if expsum == self.targets:
                coeffs += reduce(mul, [m.nchoosekm(s) for m, s in zip(self.multinoms, seq)])

        return coeffs

    def __str__(self):
        #First we need to sort the multinomials by their exponent.
        sortedmns = sorted(self.multinoms, key=(lambda m: (m.exponent,m.power)), reverse=True)
        return str(self.coefficient) + ''.join([str(mn) for mn in sortedmns])            

class Multinomial(object):
    """Represents a multinomial expansion.

    Attributes:
        power (list): The power on each of the *unexpanded* variables in the multinomial;
          of the form (x^2+y^2+z^2) => 2.

        coef (list): The coefficient on each of the variables in the multinomial; e.g.
          (x^2 + a y^2) => a. For more than two variables, each color that *has* arrowing
          gets the coefficient while the others get 1.

        arrowings (list): The exponents of the arrows terms.
        exponent (int): The exponent of the entire multinomial, default value is 1.

        powersum (list): The integer value that all term exponents in the multinomial should
          sum to (or be less than).
    
        possible_powers (list): The possible powers based on the exponent in the multinomial.
    """
    def __init__(self, power, coeff, arrowings, exponent=1):
        """Sets up the multinomial.
        
        Args:
            powers (list): The power on each of the *unexpanded* variables in the multinomial;
              of the form (x^2+y^2+z^2) => 2.

            coeff (list): The coefficient on each of the variables in the multinomial; e.g.
              (x^2 + a y^2) => a. For more than two variables, each color that *has* arrowing
              gets the coefficient while the others get 1.

            exponent (int, optional): The exponent of the entire multinomial, default is 1.
            arrowings (list): The exponents of the arrow terms.
        """
        self.power = power
        self.coeff = coeff
        self.arrowings = arrowings
        self.exponent = exponent
        self.powersum = power*exponent
        """Returns the integer value that all term exponents in the multinomial should
        sum to (or be less than)."""
        self.possible_powers = list(range(0,power*exponent+1, power))
        """For each variable being considered, determines the possible powers based
        on the exponent in the multinomial."""

    def __str__(self):
        #We want to print the multinomial out in a nice, readable way, similar to how
        #they are presented in Mathematica.
        contents = ' + '.join(["{}^{}".format(self.coeff if a else 1, self.power) for a in self.arrowings])
        return "({})^{}".format(contents, self.exponent)

    def normed_seq(self, seq):
        """Normalizes the specified sequence using the powers of unexpanded terms in 
        the multinomial.
        
        Args:
            seq (list): A list of exponents in an *expanded* term.

        Returns:
            norm_seq (list): The normalized lust of exponents.
        """
        return [int(ai/self.power) for ai in seq]

    def nchoosekm(self, sequence):
        """Returns the number of different ways to partition an n-element
        set into disjoint subsets of sizes k1, ..., km.

        Args:
            sequence (tuple): An un-normed tuple of form (k1, k2, k3).
        
        Returns:
            nckm (int): The value of the multinomial coefficient for the series.
        """
        prod = 1
        if not all(seq%self.power == 0 for seq in sequence):
            return 0
        else:
            from operator import mul
            normseq = self.normed_seq(sequence)
            len_seq = len(sequence)
            for i in range(len_seq):
                nsum = sum(normseq[0:i+1])
                prod *= Multinomial.nchoosek(nsum, normseq[i])

            #Add the contribution from the coefficients of the variable *inside*
            #the multionomial.
            pcoeff = 1
            for iseq, arrow in enumerate(self.arrowings):
                if arrow:
                    pcoeff *= self.coeff**normseq[iseq]

            return prod*pcoeff
        
    @staticmethod
    def nchoosek(n, k):
        """This implementation was taken from "Binomial Coefﬁcient Computation: Recursion 
        or Iteration?" by Yannis Manolopoulos, ACM SIGCSE Bulletin InRoads, Vol.34, No.4, 
        December 2002. http://delab.csd.auth.gr/papers/SBI02m.pdf It is supposed to be robust 
        against large, intermediate values and to have optimal complexity.

        Args:
            n (int): The value of n in n choose k.
            k (int): The value of k in n choose k.

        Return:
            t (int): The value of n choose k.
        """
        if not (-1 < k < n+1): #pragma: no cover
            return 0
        if k==0 and n == 0:
            return 1
        t = 1
        if k < n-k:
            for i in range(n, n-k, -1):
                t = t*i/(n-i+1)
        else:
            for i in range(n, k, -1):
                t = t*i/(n-i+1)

        return t

def group(gen): #pragma: no cover
    """Generates an entire group using the specified generators by applying generators
    to each other recursively.

    Args:
        gen (list): A list of generators as integers.

    Returns:
        groupi (list): A list of the symmetry group operators.
    """
    from operator import itemgetter as iget
    def g_apply(operations, source, groupi=None):
        """Applies the specified group operations to the source list of elements and then
        appends it to the group if it is unique.

        Args:
            operations (list): The symmetry group operation.
            source (list): Object for group to act on.
            groupi (list, optional): The growing symmetry group. Default is None.

        Returns:
            result (list): The result of the symmetry group on the object.
        """
        result = list(iget(*operations)(source))
        if groupi is not None and result not in groupi:
            groupi.append(result)
        return result

    #Make sure the group is zero-based for python.
    if not 0 in gen[0]:
        ngens = [list([e-1 for e in g]) for g in gen]
    else:
        ngens = gen

    groupi = []
    for i in ngens:
        for j in ngens: #filter(lambda k: k!=i, ngens):
            c = g_apply(i, j, groupi)
            d = g_apply(i, c, groupi)
            while d != c:
                d = g_apply(i, d, groupi)

    while True:
        group2 = []
        for i in ngens:
            for h in groupi:
                d = g_apply(i, h)
                if d not in groupi and d not in group2:
                    group2.append(d)

        groupi.extend(group2)
        if len(group2) == 0:
            break
    return(groupi)

def _group_to_cyclic(group, limit=None):
    """Determines the degeneracy of each r-cycle in the specified group operations.

    Args:
        group (list): The a list of the symmetry group operations.
        limit (list, optional): The the start and finishing indices of the desired group.
          Deafult None.

    Returns:
        result (list): The group converted to a cyclic group.
    """
    result = []
    #We allow filtering so that the unit testing can access the cyclic form of the group.
    if 0 not in group[0][0]: #pragma: no cover
        group = [[[j - 1 for j in i] for i in t] for t in group]
    
    if limit is not None: #pragma: no cover
        filtered = group[limit[0]:limit[1]]
    else:
        filtered = group

    for operation in filtered:
        #visitedp has the same # of elements as the site group operation and
        #is used to make sure each element in the array is visited as
        #we loop through in a *non-sequential* order.
        #visitedc has the same number of elements as the arrow group
        #opertaion and is used to make sure each element of that array
        #is visited as we loop through it in a "non-sequential" order.
        visitedp = [0]*len(operation[0])
        visitedc = [0]*len(operation[1])
        polynomials = {}
        
        arrow_cycle_len = []
        #We start by finding the cylce lengths of the arrow group
        #operations so that we can later see how many of them are
        #divosors of the site group operations.
        while 0 in visitedc:
            #Start with the first elemenet in the arrow group that
            #hasn't been visited yet. All cycles have length > 0.
            cursor = vindex = visitedc.index(0)
            cyc_len = 1
            visitedc[cursor] = 1
            cursor = operation[1][cursor]
            while cursor != vindex:
                visitedc[int(cursor)] = 1
                cyc_len += 1
                cursor = int(operation[1][int(cursor)])
            arrow_cycle_len.append(cyc_len)

        while 0 in visitedp:
            #Start with the first element in the site group that hasn't
            #been visited yet. The first non-trivial polynomials have
            #powers > 0.
            cursor = vindex = visitedp.index(0)
            powers = 1 
            #change the current position to having been visited; move the cursor.
            visitedp[cursor] = 1
            cursor = operation[0][cursor]
            #The power of the variables in the polynomials is equal to
            #the number of group operations separating the cursor's
            #current position from its *value* in the group operations
            #list.
            while cursor != vindex:
                visitedp[cursor] = 1
                powers += 1
                cursor = operation[0][cursor]
            #We now have everything need to construct part of the
            #polynomial. This is done by taking powers and using it to
            #construct an array of length equal to the number of
            #elements in the system each entry in the array is set to
            #be equal to powers.

            #To get the coefficients right we need to see how many of
            #the arrow cycles are divisors in length of the current
            #cycle.
            coef = sum([l for l in arrow_cycle_len if powers%l == 0])
            polykey = (powers, coef)
            if polykey not in polynomials:
                polynomials[polykey] = 1
            else:
                polynomials[polykey] += 1
        result.append(polynomials)
        
    return result

def polya(concentrations, group, arrowings=None, debug=False):
    """Uses a group and concentrations to find the number of unique arrangements as described by 
    polya.

    Args:
        concentrations (list):: A list of integers specifying how many of each coloring should
          be present in each of the enumerated lists.

        group (list): A list of group operations for permuting the colorings.
        arrowings (int, optional): The number of arrows present. Default is None.
        debug (bool, optional): True if the code is being debugged.

    Returns:
        polya (int): The polya number for the system.

    Raises:
        ValueError: if the concentrations don't sum to the size of the group operation.
    """
    import itertools
    
    temp = [str(i) for i in concentrations]
    for i in group:
        temp = [str(j) for j in i[0]]

    for i in group:
        temp = [str(j) for j in i[1]]
    
    if arrowings is None:
        arrowings = [False]*len(concentrations)
    elif isinstance(arrowings, int):
        arrowings = [False]*(len(concentrations)-arrowings) + [True]*arrowings
        
    #This is to check that the concentrations sum to the number of sites the group is
    #operating on
    if sum(concentrations) != len(group[0][0]):
        raise ValueError("The concentrations don't sum to the number of things the group is acting on!")
            
    len_g = len(group)
    for k, g in itertools.product(range(2),range(len_g)):
        if 0 not in group[g][k]:
            group[g][k] = [j-1 for j in group[g][k]]

    polyndict = {}
    #The operations in the group are used to construct the unique polynomials for each operation.
    for polynomials in _group_to_cyclic(group):
        #Construct a product of multinomials for this group operation.
        p = Product(1,concentrations)
        for exp, coeff in polynomials:
            p.multinoms.append(Multinomial(exp, coeff, arrowings, polynomials[(exp, coeff)]))

        key = str(p)
        if key not in polyndict:
            polyndict[key] = p
        else:
            polyndict[key].coefficient += 1

    if debug: #pragma: no cover
        for key in polyndict:
            print((str(polyndict[key]), " => ", polyndict[key].coeff()))
            
    pp = 0
    for p in list(polyndict.values()):
        pp += p.coeff()
        
    rad = sum([p.coeff() for p in list(polyndict.values())])
    return int(rad/float(len(group)))
