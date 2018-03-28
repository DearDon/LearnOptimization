program Newton
	real*8 x(2),t,g(2),A(2,2),Gni(2,2),b(2),f,p(2),ans,aro(2)
	data x /0,0/
	data A /2,-1,-1,2/
	data Gni /0.6666667,0.3333334,0.3333334,0.666666667/
	data b /-10,-4/
	open(13,file="Newton")
	write(13,"(5A13)")"f","x1","x2","g1","g2"
	do
	g=matmul(A,x)+b
	f=60-10*x(1)-4*x(2)+x(1)*x(1)+x(2)*x(2)-x(1)*x(2)
	write(13,"(5f13.6)")f,x,g
	call onedimenmul(g,g,2,ans)
	if(ans<0.01) exit
	p=matmul(-Gni,g)
	x=x+p
	end do
	close(13)
end program Newton

	subroutine onedimenmul(m1,m2,n,ans)
	integer n
	real*8 m1(n),m2(n),ans
	ans=0
	do i=1,n
		ans=m1(i)*m2(i)+ans
	end do
	end subroutine onedimenmul
