!
      subroutine myspline(n,x,y,nn,xn,yn)
!c---------spline fit to derive the yy value at point xx
!c  Inputs:
!c     n:    the length of x and y
!c     x(n): the x values which x(1) < x(2) ... < x(n)
!c     y(n): the y value which correspondent to x(n)
!c     nn:  the length of vector xx and yy
!c     xn:  the x value at which y value is wanted
!c
!c  Outputs:
!c     yn: the wanted y value from the fitting
!c
!c  Internal variables:
!c     yp1: the derivative of y over x at x(1), for natural bc, yp1=1.e31
!c     ypn: the derivative of y over x at x(n), for natural bc, ypn=1.e31
!c     y2(n): the second derivatives
      IMPLICIT NONE
      integer :: n,nn,i
      INTEGER, PARAMETER :: ny2=5000
      real(kind=8) :: x(n),y(n),xn(nn),yn(nn),y2(ny2),xx,yy,yp1,ypn
!c--------the sorting which makes sure x(1)<x(2)<...<x(n)-------
        call sort2(n,x,y) 
!c--------start spline------------
        yp1=1.e31
        ypn=1.e31
          call spline(x,y,n,yp1,ypn,y2)
          do i=1,nn
          xx = xn(i)
          IF( xx > x(1) .and. xx < x(n) ) THEN
          call splint(x,y,y2,n,xx,yy)
          yn(i) = yy
          ELSE IF( xx <= x(1) ) THEN
            yn(i) = y(1)
          ELSE
            yn(i) = y(n)
          END IF
         enddo
         return
         end 
!
      SUBROUTINE spline(x,y,n,yp1,ypn,y2) 
      IMPLICIT NONE
      INTEGER n
      INTEGER, PARAMETER :: NMAX=5000
      INTEGER, PARAMETER :: ny2=5000
      real(kind=8) :: yp1,ypn,x(n),y(n),y2(ny2) 
      INTEGER i,k 
      real(kind=8) :: p,qn,sig,un,u(NMAX) 
      if (yp1.gt..99e30) then 
        y2(1)=0. 
        u(1)=0. 
      else 
        y2(1)=-0.5 
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1) 
      endif 
      do 11 i=2,n-1 
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1)) 
        p=sig*y2(i-1)+2. 
        y2(i)=(sig-1.)/p 
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)      &
    -x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*    & 
     u(i-1))/p 
11    continue 
      if (ypn.gt..99e30) then 
        qn=0. 
        un=0. 
      else 
        qn=0.5 
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1))) 
      endif 
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.) 
      do 12 k=n-1,1,-1 
        y2(k)=y2(k)*y2(k+1)+u(k) 
12    continue 
      return 
      END 
!
      SUBROUTINE splint(xa,ya,y2a,n,x,y) 
      IMPLICIT NONE
      INTEGER n 
      real(kind=8) :: x,y,xa(n),y2a(n),ya(n) 
      INTEGER k,khi,klo 
      REAL a,b,h 
      klo=1 
      khi=n 
1     if (khi-klo.gt.1) then 
        k=(khi+klo)/2 
        if(xa(k).gt.x)then 
          khi=k 
        else 
          klo=k 
        endif 
      goto 1 
      endif 
      h=xa(khi)-xa(klo) 
      if (h.eq.0.) then
        print *, 'bad xa input in splint' 
        stop
      endif
      a=(xa(khi)-x)/h 
      b=(x-xa(klo))/h 
       y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+ &
        (b**3-b)*y2a(khi))*(h*h)/6. 
      return 
      END 
!C  (C) Copr. 1986-92 Numerical Recipes Software 71a. 
      SUBROUTINE sort2(n,arr,brr) 
      IMPLICIT NONE
      INTEGER n,M,NSTACK 
      real(kind=8) :: arr(n),brr(n) 
      PARAMETER (M=7,NSTACK=50) 
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK) 
      REAL a,b,temp 
      jstack=0 
      l=1 
      ir=n 
1     if(ir-l.lt.M)then 
        do 12 j=l+1,ir 
          a=arr(j) 
          b=brr(j) 
          do 11 i=j-1,1,-1 
            if(arr(i).le.a)goto 2 
            arr(i+1)=arr(i) 
            brr(i+1)=brr(i) 
11        continue 
          i=0 
2         arr(i+1)=a 
          brr(i+1)=b 
12      continue 
        if(jstack.eq.0)return 
        ir=istack(jstack) 
        l=istack(jstack-1) 
        jstack=jstack-2 
      else 
        k=(l+ir)/2 
        temp=arr(k) 
        arr(k)=arr(l+1) 
        arr(l+1)=temp 
        temp=brr(k) 
        brr(k)=brr(l+1) 
        brr(l+1)=temp 
        if(arr(l+1).gt.arr(ir))then 
          temp=arr(l+1) 
          arr(l+1)=arr(ir) 
          arr(ir)=temp 
          temp=brr(l+1) 
          brr(l+1)=brr(ir) 
          brr(ir)=temp 
        endif 
        if(arr(l).gt.arr(ir))then 
          temp=arr(l) 
          arr(l)=arr(ir) 
          arr(ir)=temp 
          temp=brr(l) 
          brr(l)=brr(ir) 
          brr(ir)=temp 
        endif 
        if(arr(l+1).gt.arr(l))then 
          temp=arr(l+1) 
          arr(l+1)=arr(l) 
          arr(l)=temp 
          temp=brr(l+1) 
          brr(l+1)=brr(l) 
          brr(l)=temp 
        endif 
        i=l+1 
        j=ir 
        a=arr(l) 
        b=brr(l) 
3       continue 
          i=i+1 
        if(arr(i).lt.a)goto 3 
4       continue 
          j=j-1 
        if(arr(j).gt.a)goto 4 
        if(j.lt.i)goto 5 
        temp=arr(i) 
        arr(i)=arr(j) 
        arr(j)=temp 
        temp=brr(i) 
        brr(i)=brr(j) 
        brr(j)=temp 
        goto 3 
5       arr(l)=arr(j) 
        arr(j)=a 
        brr(l)=brr(j) 
        brr(j)=b 
        jstack=jstack+2 
        if(jstack.gt.NSTACK) then
          print *, 'NSTACK too small in sort2' 
          STOP
        endif
        if(ir-i+1.ge.j-l)then 
          istack(jstack)=ir 
          istack(jstack-1)=i 
          ir=j-1 
        else 
          istack(jstack)=j-1 
          istack(jstack-1)=l 
          l=i 
        endif 
      endif 
      goto 1 
      END 
!C  (C) Copr. 1986-92 Numerical Recipes Software 71a. 
