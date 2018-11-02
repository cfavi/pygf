import math

class GF2n:
    '''
    Model GF(2**n) elements.

    Objects of this class act as elements of GF(2**n). 

    >>> GF2n(0x11b, b'\\x11\\xb0\\x00\\x01') == GF2n(0x11b, 1)
    True
    '''
    def __init__(self, poly, val=0):
        if isinstance(poly, (bytes, bytearray)):
            self.p = int.from_bytes(poly, byteorder='big')
        else:
            self.p = poly
        self.n = self.p.bit_length() - 1 # len(bin(poly))-2-1

        if isinstance(val, (bytes, bytearray)):
            self.val = int.from_bytes(val, byteorder='big')
        else:
            self.val = val
        self._reduce_modp()
        self.log = False

    def __repr__(self):
        return 'GF2n(0x{:x}, 0x{:x})'.format(self.p, self.val)

    def __str__(self):
        ''' Polynomial representation
        >>> print(GF2n(0x11b, 3))
        x + 1 (0x3)
        >>> print(GF2n(0x11b, 0x1b))
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
        '''Redefine equality. '''
        if isinstance(b, GF2n):
            return (self.val == b.val) and (self.p == b.p)
        else:
            return self == GF2n(self.p, b)

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
    
    def __sub__(self, b):
        ''' Addition in Galois Field is XOR
        >>> a = GF2n(0x11b, 3) 
        >>> b = GF2n(0x11b, 1)
        >>> b - a == GF2n(0x11b, 2)
        True
        '''
        return self.__add__(b)


    def __mul__(self, b):
        ''' GF Multiplication double and add implementation
        >>> a = GF2n(0x11b, 3)
        >>> a * a == GF2n(0x11b, 5)
        True
        '''
        if b == 1:
            return GF2n(self.p, self.val)
        elif b == 2:
            t = self.val * 2
            if t >= 2**self.n:  # quick reduction
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
        ''' Exponentiation 
        >>> GF2n(0x11b, 3) ** 2 == GF2n(0x11b, 5)
        True
        '''
        if not isinstance(other, int):
            return NotImplemented

        #bit scanning implementation
        # r = GF2n(self.p, 1)
        # p = self * 1
        # while other:
        #     if other & 1:
        #         r *= p
        #     p = p * p
        #     other >>= 1
        # return r

        #montgomery ladder implementation
        x0 = GF2n(self.p, 1)
        x1 = GF2n(self.p, self.val)
        for kj in bin(other)[2:]:
            if kj == '0':
                x1 = x0 * x1
                x0 = x0 * x0
            else:
                x0 = x0 * x1
                x1 = x1 * x1
            if self.log: print("x0={}\tx1={}".format(x0,x1))
        return x0
    
    def __bytes__(self):
        '''Create a serialized bytes string of the value (and not the poly)

        The byte string is MSBf coded on math.ceil(self.n/8) bytes
        >>> bytes(GF2n(0x11b, 5)) 
        b'\\x05'
        >>> bytes(GF2n(0x1277, 5))
        b'\\x00\\x05'
        '''
        return self.val.to_bytes(math.ceil(self.n/8), byteorder='big')
        
    def _reduce_modp(self):
        '''Reduce self.val by self.p so that the new value is < 2**self.n

        >>> GF2n(0x11b, 0x11b) == GF2n(0x11b, 0)
        True
        >>> GF2n(0x11b, 0x100) == GF2n(0x11b, 0x1b)
        True
        >>> GF2n(0x11b, 0x11b027) == GF2n(0x11b, 0x27)
        True
        '''
        if self.val < 2**self.n:
            return  # Do nothing

        while self.val >= 2**self.n:
            valn = self.val.bit_length() - 1  # len(bin(self.val))-2-1
            self.val ^= self.p << (valn - self.n)
        


if __name__ == '__main__':
    import doctest
    doctest.testmod()

