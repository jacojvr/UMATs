      subroutine shearmod(eg,enu,temp,props,nprops)
      implicit real*8(a-h,o-z)
      dimension props(nprops)
      emu0=props(1)
      ed0=props(2)
      et0=props(3)
      enu=min(dabs(props(4)),0.499d0) 
      if(temp.gt.et0) then
        eg = emu0 - ed0/(dexp(et0/temp)-1.d0)
      else
        eg = emu0
      endif
      return
      end