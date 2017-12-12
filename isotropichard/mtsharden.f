      subroutine fisotropic(sy,dsy,deqpl,dtime,temp,
     & statev,tempstatev,nstatv,props,nprops)
      implicit real*8(a-h,o-z)
      dimension statev(nstatv),tempstatev(nstatv),props(nprops)
c
      character*8 cmname
      parameter(zero=0.d0,one=1.d0,two=2.d0,newtonmax=10,toler=1.0d-8,
     & ratemin=1.d-10)
c
c     material properties
      emu0=props(1)
      ed0=props(2)
      et0=props(3)
c     props(4) is used for poisson's ratio
      sa=props(5)
      se0s=props(6)
      si=props(7)
      a0e=props(8)
      a0es=props(9)
      a0i=props(10)
      hard0=props(11)
      rate0=props(12)
      rate0es=props(13)
      qe=props(14)
      pe=props(15)
      qi=props(16)
      pi=props(17)
      pow=props(18)
      cm=props(19)
c
      cmname = 'MTS_POWR'
      if(cm.gt.zero)then
        cmname = 'MTS_TANH'
      endif
c
c     check 
      rate = deqpl/dtime
      if(rate.lt.ratemin)then
        rate = ratemin
      endif
c
c     harden in statev
      sep = statev(1)
      sec = sep
c
      emu = emu0 - ed0/(dexp(et0/temp)-one)
      emuSF = emu/emu0
c  
      sfi0 = temp/(a0i*emu)
      sfli = dlog(rate0/rate)*sfi0
      sfi = (one-sfli**(one/qi))**(one/pi)
c
      sfe0 = temp/(a0e*emu)
      sfle = dlog(rate0/rate)*sfe0
      sfe = (one-sfle**(one/qe))**(one/pe)
c    
      c0 = -temp/(emu*a0es)
      sat = se0s*(rate/rate0es)**c0
c    
      if(cmname(5:8).eq.'TANH')then
        hard = hard0*(one-(dtanh(pow*sec/sat)/dtanh(pow)))
        dhard = deqpl*hard0*pow/
     &   (dtanh(pow)*sat*dcosh(pow*sec/sat)**two)
      else
        hard = hard0*(one-(sec/sat))**pow
        dhard = (deqpl*hard0*pow*(one-(sec/sat))**(pow-one))/sat
      endif
      f = sec - sep - deqpl*hard
c
      if(dabs(f).gt.toler)then
      kount=1
      do while((dabs(f).gt.toler).and.(kount.lt.newtonmax))
        kount=kount+1
        sec = sec - f/(one+dhard)
        if(cmname(5:8).eq.'TANH')then
          hard = hard0*(one-(dtanh(pow*sec/sat)/dtanh(pow)))
          dhard = deqpl*hard0*pow/
     &            (dtanh(pow)*sat*dcosh(pow*sec/sat)**two)
        else
          hard = hard0*(one-(sec/sat))**pow
          dhard = (deqpl*hard0*pow*(one-(sec/sat))**(pow-one))/sat
        endif  
        f = sec - sep - deqpl*hard
      enddo
      endif
c
c     Yield stress
      sy = sa + emusf*(sfi*si+sfe*sec)
c
c     Partial gradient componenets
c     d(sec)/d(epl)
      dsecdeqpl = hard/(one+dhard)  
c     d(sfi)/d(epl)
      dsfideqpl = sfi0*(one-sfli**(one/qi))**(one/pi-one)*
     &          sfli**(one/qi-one)/(pi*qi*rate)
c     d(sfe)/d(epl)
      dsfedeqpl = sfe0*(one-sfle**(one/qe))**(one/pe-one)*
     &          sfle**(one/qe-one)/(pe*qe*rate)
c     d(sec)/d(rate)
      dsecdrt = c0*dhard*se0s*sec*(rate/rate0es)**(c0-one)/
     &          ((dhard+one)*rate0es*sat)
c     Total
c     d(yield)/d(epl)
      dsy =  emuSF*(sfe*(dsecdeqpl+dsecdrt/dtime)+dsfedeqpl*sec/dtime+
     &       dsfideqpl*si/dtime)
c
      tempstatev(1) = sec
      return
      end