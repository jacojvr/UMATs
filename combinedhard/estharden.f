      subroutine fcombined(sy,dsy,sb,dsb,deqpl,dtime,temp,
     & statev,tempstatev,nstatv,props,nprops)
c
      implicit real*8(a-h,o-z)
      dimension statev(nstatv),tempstatev(nstatv),props(nprops)
c
      parameter(zero=0.d0,half=0.5d0,one=1.d0,two=2.d0,newtonmax=10,
     & toler=1.0d-8,ratemin=1.d-10)
c
c     material properties
      emu0=props(1)
      ed0=props(2)
      et0=props(3)
c     props(4) is used for poisson's ratio
c     stress values
      sa=props(5)
      s0=props(6) 
c     scaling function
      a0e=props(7)
      pe=props(8)
      qe=props(9)
      rate0=props(10)
c     evolve misorient
      cld=props(11)
c     geometric disl
      cg=props(12)
      rg=props(13)
c     accumulation
      c1=props(14)
c     dynamic recovery
      c20=props(15)
      a02=props(16)
      rate02=props(17)
c     thermal recovery
      c30=props(18)
      a03=props(19)
      r3=props(20)
c     recoverable
      c4=props(21)
      c50=props(22)
      a05=props(23)
c     back stress
      c6=props(24)
      c7=props(25)
c
c     check 
      rate = deqpl/dtime
      if(rate.lt.ratemin)then
        rate = ratemin
      endif
c
c     statevariables
      x1p=statev(1) !misorientation
      x2p=statev(2) !dislocation ratio
      x3p=statev(3) !recoverable disl.
c     previous back stress:
      sbp=statev(4)
c
      emu=emu0-ed0/(dexp(et0/temp)-one)
      emusf=emu/emu0
c
c     scaling function
      sfe0=temp/(a0e*emu)
      sfle=dlog(rate0/rate)*sfe0
      sfe=(one-sfle**(one/qe))**(one/pe)
c
c     dynamic recovery
      c2t=-temp/(emu*a02)
      c2=c20*dexp(c2t*dlog(rate/rate0))
      dc2drt=c2*c2t
c
c     thermal recovery
      c3=c30*dexp(-a03/temp)*dtime
      c5=c50*dexp(-a05/temp)*dtime
c    
c     update lattice incompatibility
      x1=x1p+cld*deqpl
c     initial guess for dislocation dens.
      x2=x2p
      sqx2=dsqrt(x2)
c
      f=-deqpl*(cg*x1**rg+c1*sqx2-c2*x2)+c3*x2**r3
      ddf=one-deqpl*(half*c1/sqx2-c2)
      if(dabs(f).gt.toler)then
      kount=0
      do while((dabs(f).gt.toler).and.(kount.lt.newtonmax))
      kount=kount+1
      x2=dabs(x2-f/ddf)
      sqx2=dsqrt(x2)
      f=x2-x2p-deqpl*(cg*x1**rg+c1*sqx2-c2*x2)+c3*x2**r3
      ddf=one-deqpl*(half*c1/sqx2-c2)
      enddo
      endif
      x3=(x3p+deqpl*c4*sqx2)/(one+deqpl*c2+c5)
c
c     equivalent threshold and yield stress
      sec=s0*sqx2
      sy=sa+emusf*sfe*sec
c     back stress
      sbnum=sbp+deqpl*c6*sqx2
      sbden=one+deqpl*c7
      sb=sbnum/sbden
c
c     dislocation sensitivity (EQ. 5.44)
      dx2de=(cg*x1**rg+c1*sqx2-c2*x2-dc2drt*x2+
     & deqpl*rg*cg*cld*(x1**(rg-one)))/ddf
c     isotropic senistivity
      dsec=half*s0*dx2de/sqx2
      dsfe=sfe0*(one-sfle**(one/qe))**(one/pe-one)*
     & sfle**(one/qe-one)/(pe*qe*rate*dtime)
      dsy=emusf*(sfe*dsec+dsfe*sec)
c     back stress sensitivity
      dsb=(c6*sqx2 + half*deqpl*c6*dx2de/sqx2)/sbden-
     & c7*sbnum/(sbden**two)
c     
c     update temporary state variable array
      tempstatev(1)=x1
      tempstatev(2)=x2
      tempstatev(3)=x3
      tempstatev(4)=sb
c
      return
      end