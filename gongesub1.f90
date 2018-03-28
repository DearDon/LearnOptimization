module gongesub1
	contains 
!The following subroutine is the to get x from Ax=b by Conjugate Gradient
!N is the N_max ,this method is come from "矩阵计算的理论与方法"北大徐树方p154 
!with a little change
	subroutine cg(A,b,x,N)
	!variable	meaning
	!A		a matrix,from main program,is coefficient matrix of Ax=b
	!b		a vector,from main program,is righthand of Ax=b
	!x		a vector,the answer of Ax=b,is what we need,our goal
	!r		a vector,minus grads of 0.5xAx-bx at x point,says b-Ax
	!p		a vector,the direction of iteration better than r
	!w		a vector,value is A*p,is useful to simplify the process
	!q0		a number,value is r0*r0,is standard of loop times
	!q1		a number,value is rk-1*rk-1
	!q2		a number, value is rk*rk
	!ba,ar		a number,named by their pronounciation
	!e		a number,standard of loop times,input by client
	!test		a number,value is matmul(r,w)
	!pw		a number,value is matmul(p,w)
	!i		a number,count variable
	!N		a number,the degree of A
		real*8 A(N,N)
		real*8 b(N),x(N),r(N),p(N),w(N)
		!real*8 A(2,2),b(2),x(2),r(2),p(2),w(2)
		real*8 q0,q1,q2,ba,ar,e,test,pw
		integer i,N
	!	write(*,*)"you want the x_error less than"
	!	read(*,*)e
	!	write(*,*)"you want x0 to be?"
	!	read(*,*)x
		e=0.01d0
		r=b-matmul(a,x)
		call onedimenmul(r,r,N,q0)
		q2=q0
		p=r

	!	w=matmul(A,p)
	!	call onedimenmul(r,w,2,test)
	!	ar=q2/test
	!	x=x+ar*p
	!	r=r-ar*w
	!	q1=q2;call onedimenmul(r,r,2,q2)
	!	!r=r-a*w
		i=1
	write(*,"(5f13.6,i13)")2*x(1)**2+x(2)**2-x(1)*x(2),x,r,i
	
		do while(q2>=e)
		q1=q2
	!	ba=q2/q1
	!	p=r+ba*p
		
		w=matmul(A,p)
		call onedimenmul(p,w,N,pw)!pw is p*w
		ar=q1/pw
		x=x+ar*p
		r=r-ar*w
		call onedimenmul(r,r,N,q2)
		!r=r-a*w
		ba=q2/q1
		i=i+1
		p=r+ba*p
	write(*,"(5f13.6,i13)")2*x(1)**2+x(2)**2-x(1)*x(2),x,r,i
	
		end do
!	write(*,*)"x",x
!	write(*,*)"i",i
	end subroutine cg

	!This subroutine is to solve one dimention's multiplication
	subroutine onedimenmul(m1,m2,n,ans)
	integer n
	real*8 m1(n),m2(n),ans
	ans=0
	do i=1,n
		ans=m1(i)*m2(i)+ans
	end do
	
	end subroutine onedimenmul
end module gongesub1

