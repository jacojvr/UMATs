      subroutine shearmod(eg,enu,temp,props,nprops)
      implicit real*8(a-h,o-z)
      dimension props(nprops)
      eg = props(1)
      enu=min(dabs(props(2)),0.4999d0)
      return
      end