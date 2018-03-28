program DFP
	real*8 x(2),t,g(2),A(2,2),f,p(2),b(2),ans,aro(2),s(2),y(2),H(2,2),Hf(2,2)
	data x /8,9/
	data A /8,0,0,2/
	data b /40,12/
	data H /1,0,0,1/
	open(13,file="5.6.data")
	write(13,"(5A13)")"f","x1","x2","g1","g2"
	do
	g=matmul(A,x)-b
	f=4.0*(x(1)-5)**2+(x(2)-6)**2
	write(13,"(5f13.6)")f,x,g
	call onedimenmul(g,g,2,ans)
	if(ans<0.01) exit
	p=matmul(-H,g)
	call onedimenmul(p,g,2,t)
	aro=matmul(a,p)
	call onedimenmul(p,aro,2,ans)
	t=t/ans
	y=matmul(a,x-t*p)-matmul(a,x)
	s=t*p
	x=x-t*p
	do i=1,2
	do j=1,2
	hf(i,j)=s(i)*s(j)
	end do
	end do
	call onedimenmul(s,y,2,ans)
	h=h+(1/ans)*hf
	aro=matmul(h,y)
	do i=1,2
	do j=1,2
	hf(i,j)=aro(i)*aro(j)
	end do
	end do
	call onedimenmul(y,aro,2,ans)
	h=h-(1/ans)*hf
	
	end do
	close(13)

end program DFP

	subroutine onedimenmul(m1,m2,n,ans)
	integer n
	real*8 m1(n),m2(n),ans
	ans=0
	do i=1,n
		ans=m1(i)*m2(i)+ans
	end do
	end subroutine onedimenmul
