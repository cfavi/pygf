

class GF2n:
    def __init__(self, poly, val=0):
        self.n = len(bin(poly))-2-1
        self.p = poly
        self.val = val #should be reduced by p

    def __repr__(self):
        return 'GF2n(0x{:x}, 0x{:x})'.format(self.p, self.val)

    def __str__(self):
        ''' Polynomial representation
        >>> print GF2n(0x11b, 3)
        x + 1 (0x3)
        >>> print GF2n(0x11b, 0x1b)
        x^4 + x^3 + x + 1 (0x1b)
        '''
        s = []
        for i in reversed([x[0] for x in enumerate(reversed(bin(self.val)[2:])) if x[1]=='1']):
            if (i==0):
                s.append('1')
            elif (i==1):
                s.append('x')
            else:
                s.append('x^{}'.format(i)) 
        if len(s) == 0:
            s.append('0')
        __str = ' + '.join(s)
        __str += ' (0x{:x})'.format(self.val)
        return __str

    def __eq__(self, b):
        ''' redefine equality '''
        if isinstance(b, GF2n):
            return (self.val == b.val) and (self.p == b.p)
        else:
            return self.val == b

    def __add__(self, b):
        ''' Addition in Galois Field is XOR
        >>> a = GF2n(0x11b, 3) 
        >>> b = GF2n(0x11b, 1)
        >>> a + b == GF2n(0x11b, 2)
        True
        '''
        if not isinstance(b, GF2n):
            b = GF2n(self.p, b)
        return GF2n(self.p, self.val ^ b.val)


    def __mul__(self, b):
        ''' GF Multiplication dumb implementation
        >>> a = GF2n(0x11b, 3)
        >>> a * a == GF2n(0x11b, 5)
        True
        '''
        if b == 1:
            return GF2n(self.p, self.val)
        elif b == 2:
            t = self.val * 2
            if t >= 2**self.n:
                t ^= self.p
                
            return GF2n(self.p, t)

        if not isinstance(b, GF2n):
            b = GF2n(self.p, b)
        
        m = b.val
        p = self * 1
        r = GF2n(self.p, 0)
        while m:
            if m & 1:
                r += p
            p = p * 2
            m = m >> 1
        
        return r

    def __rmul__(self, a):
        ''' right mul to be able to do 2 * GF2n(0x11b, 0x34)'''
        return self * a
        

    def __pow__(self, other):
        ''' exponentiation 
        >>> GF2n(0x11b, 3) ** 2 == GF2n(0x11b, 5)
        True
        '''
        if not isinstance(other, int):
            return NotImplemented
        
        r = GF2n(self.p, 1)
        p = self * 1
        while other:
            if other & 1:
                r *= p
            p = p * p
            other >>= 1

        return r