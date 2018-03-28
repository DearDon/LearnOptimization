program Newtonxiu
	real*8 x(2),t,g(2),A(2,2),Gni(2,2),b(2),f,p(2),ans,aro(2)
	data x /0,0/
	data A /8,0,0,4/
	data Gni /0.125,0.0,0.0,0.25/
	data b /9,-3/
	open(13,file="Newtonxiu")
	write(13,"(5A13)")"f","x1","x2","g1","g2"
	do
	g=matmul(A,x)+b
	f=10+x(1)+x(2)+4*(x(1)+1)**2+2*(x(2)-1)**2
	write(13,"(5f13.6)")f,x,g
	call onedimenmul(g,g,2,ans)
	if(ans<0.01) exit
	p=matmul(-Gni,g)
	aro=matmul(p,A)
	t=ans
	call onedimenmul(aro,g,2,ans)
	t=t/ans
	x=x-t*p
	end do
	close(13)
end program Newtonxiu

	subroutine onedimenmul(m1,m2,n,ans)
	integer n
	real*8 m1(n),m2(n),ans
	ans=0
	do i=1,n
		ans=m1(i)*m2(i)+ans
	end do
	end subroutine onedimenmul
