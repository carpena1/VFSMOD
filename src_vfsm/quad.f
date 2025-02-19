       SUBROUTINE QUAD
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     THE SUBROUTINE QUAD DEFINES THE VALUES OF THE PARAMETERS                C
C     REQUIRED FOR THE NUMBERICAL INTEGRATION OF ELEMENT MATRICES             C
C     AND VECTORS.  THESE DATA ARE PROBLEM INDEPENT AND ARE GIVEN             C
C     OVER THE INTERVAL [-1,1].                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       COMMON/CINT/XI(20,20),W(20,20)


C-------Gaussian Quadrature of order 1-----------

       XI(1,1) = 0.D0
       W(1,1)  = 2.D0

C-------Gaussian Quadrature of order 2-----------

       XI(1,2) = -1.D0/DSQRT(3.D0)
       XI(2,2) = -XI(1,2)
       W(1,2)  = 1.D0
       W(2,2)  = W(1,2)

C-------Gaussian Quadrature of order 3-----------

       XI(1,3) = -DSQRT(3.D0/5.D0)
       XI(2,3) = 0.D0
       XI(3,3) = -XI(1,3)
       W(1,3)  = 5.D0/9.D0
       W(2,3)  = 8.D0/9.D0
       W(3,3)  = W(1,3)

C-------Gaussian Quadrature of order 4-----------

       XI(1,4) = -0.8611363116D0
       XI(2,4) = -0.3399810436D0
       XI(3,4) = -XI(2,4)
       XI(4,4) = -XI(1,4)
       W(1,4)  = 0.3478548451D0
       W(2,4)  = 0.6521451549D0
       W(3,4)  = W(2,4)
       W(4,4)  = W(1,4)

C------Gaussian quadrature of order 5------------

       XI(1,5) = -0.906179845938664D0
       XI(2,5) = -0.538469310105683D0
       XI(3,5) = 0.D0
       XI(4,5) = -XI(2,5)
       XI(5,5) = -XI(1,5)
       W(1,5) =  0.236926885056189D0
       W(2,5) =  0.478628670499366D0
       W(3,5) =  0.568888888888889D0
       W(4,5) =  W(2,5)
       W(5,5) =  W(1,5)

C------Gaussian quadrature of order 6------------
       XI(1,6) = -.932469514203252D0
       XI(2,6) = -.661209386466265D0
       XI(3,6) = -.238619186083197D0
       XI(4,6) = -XI(3,6)
       XI(5,6) = -XI(2,6)
       XI(6,6) = -XI(1,6)
       W(1,6) = 0.171324492379170D0
       W(2,6) = 0.360761573048139D0
       W(3,6) = 0.467913934572691D0
       W(4,6) =W(3,6)
       W(5,6) =W(2,6)
       W(6,6) =W(1,6)

C------Gaussian quadrature of order 7------------
       XI(1,7) = -0.949107912342759D0
       XI(2,7) = -0.741531185599394D0
       XI(3,7) = -0.405845151377397D0
       XI(4,7) = 0.D0
       XI(5,7) = -XI(3,7)
       XI(6,7) = -XI(2,7)
       XI(7,7) = -XI(1,7)
       W(1,7) = 0.129484966168870D0
       W(2,7) = 0.279705391489277D0
       W(3,7) = 0.381830050505119D0
       W(4,7) = 0.417959183673469D0
       W(5,7) =W(3,7)
       W(6,7) =W(2,7)
       W(7,7) =W(1,7)

C------Gaussian quadrature of order 8------------
       XI(1,8) = -0.960289856497536D0
       XI(2,8) = -0.796666477413627D0
       XI(3,8) = -0.525532409916329D0
       XI(4,8) = -0.183434642495650D0
       XI(5,8) = -XI(4,8)
       XI(6,8) = -XI(3,8)
       XI(7,8) = -XI(2,8)
       XI(8,8) = -XI(1,8)
       W(1,8) = 0.101228536290376D0
       W(2,8) = 0.222381034453374D0
       W(3,8) = 0.313706645877887D0
       W(4,8) = 0.362683783378362D0
       W(5,8) =W(4,8)
       W(6,8) =W(3,8)
       W(7,8) =W(2,8)
       W(8,8) =W(1,8)

C------Gaussian quadrature of order 9------------
       XI(1,9) = -0.968160239507626D0
       XI(2,9) = -0.836031107326636D0
       XI(3,9) = -0.613371432700590D0
       XI(4,9) = -0.324253423403809D0
       XI(5,9) = 0.D0
       XI(6,9) = -XI(4,9)
       XI(7,9) = -XI(3,9)
       XI(8,9) = -XI(2,9)
       XI(9,9) = -XI(1,9)
       W(1,9) = 0.081274388361574D0
       W(2,9) = 0.180648160694857D0
       W(3,9) = 0.260610696402935D0
       W(4,9) = 0.312347077040003D0
       W(5,9) = 0.330239355001260D0
       W(6,9) =W(4,9)
       W(7,9) =W(3,9)
       W(8,9) =W(2,9)
       W(9,9) =W(1,9)

C------Gaussian quadrature of order 10------------
       XI(1,10)  = -0.973906528517172D0
       XI(2,10)  = -0.865063366688985D0
       XI(3,10)  = -0.679409568299024D0
       XI(4,10)  = -0.433395394129247D0
       XI(5,10)  = -0.148874338981631D0
       XI(6,10)  = -XI(5,10)
       XI(7,10)  = -XI(4,10)
       XI(8,10)  = -XI(3,10)
       XI(9,10)  = -XI(2,10)
       XI(10,10) = -XI(1,10)
       W(1,10)  = 0.066671344308688D0
       W(2,10)  = 0.149451349150581D0
       W(3,10)  = 0.219086362515982D0
       W(4,10)  = 0.269266719309996D0
       W(5,10)  = 0.295524224714753D0
       W(6,10)  = W(5,10)
       W(7,10)  = W(4,10)
       W(8,10)  = W(3,10)
       W(9,10)  = W(2,10)
       W(10,10) = W(1,10)

C------Gaussian quadrature of order 20------------
       XI(1,20)  = -.993128599185094924776D0
       XI(2,20)  = -.963971927277913791287D0
       XI(3,20)  = -.912234428251325905857D0
       XI(4,20)  = -.839116971822218823420D0
       XI(5,20)  = -.746331906460150792634D0
       XI(6,20)  = -.636053680726515025467D0
       XI(7,20)  = -.510867001950827097985D0
       XI(8,20)  = -.373706088715419560662D0
       XI(9,20)  = -.227785851141645078076D0
       XI(10,20) = -.0765265211334973337513D0
       XI(11,20)  = -XI(10,20)
       XI(12,20)  = -XI(9,20)
       XI(13,20)  = -XI(8,20)
       XI(14,20)  = -XI(7,20)
       XI(15,20)  = -XI(6,20)
       XI(16,20)  = -XI(5,20)
       XI(17,20)  = -XI(4,20)
       XI(18,20)  = -XI(3,20)
       XI(19,20)  = -XI(2,20)
       XI(20,20) = -XI(1,20)
       W(1,20)  = 0.0176140071391521183115D0
       W(2,20)  = 0.0406014298003869413320D0
       W(3,20)  = 0.0626720483341090635663D0
       W(4,20)  = 0.0832767415767047487264D0
       W(5,20)  = .101930119817240435039D0
       W(6,20)  = .118194531961518417310D0
       W(7,20)  = .131688638449176626902D0
       W(8,20)  = .142096109318382051326D0
       W(9,20)  = .149172986472603746785D0
       W(10,20) = .152753387130725850699D0
       W(11,20)  = W(10,20)
       W(12,20)  = W(9,20)
       W(13,20)  = W(8,20)
       W(14,20)  = W(7,20)
       W(15,20)  = W(6,20)
       W(16,20)  = W(5,20)
       W(17,20)  = W(4,20)
       W(18,20)  = W(3,20)
       W(19,20)  = W(2,20)
       W(20,20) = W(1,20)

       RETURN
       END
