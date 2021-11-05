module MAT
    
    use IO
    
    implicit none
    contains

    subroutine jacobi(A,N,SD,V)

! Subrotina Jacobi
!
! Retirada do Numerical Recipes - Pag-460
!  
! Diagonaliza a matriz de Overlap (S) que nesta subrotina e o (A)
!
! Computa todos os autovalores e autovetores de uma matriz simetrica real A,
! que de tamanho NxN. No output, os elementos de A acima da diagonal sao desprezados.
! A matriz D retorna os autovalores de A em seus primeiros N elementos. V e uma matriz com 
! a mesma dimensao logica e fisica de A, na qual as colunas contem , no output, os autovetores normalizados de 
! A. 
! Nrot retorna o numero de rotacoes de Jacobi que foram requeridos.
!
! A => Matriz a ser diagonalizada.(NxN)
! D => Autovalores de A.
! V => Autovetores de A.
! SD => Matriz diagonalizada composta dos autovalores de D.
      
!   Variables            
      integer n,nrot,nmax
      real*8 a(n,n) !in
      real*8 d(n),v(n,n),sd(n,n) !out
      parameter (nmax=500)
      integer i,ip,iq,j
      real*8 c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax)

!     Inicializa  a matriz identidade

        do ip=1,n
           do iq=1,n
              v(ip,iq)=0
           enddo
           v(ip,ip)=1
        enddo
        
!     Inicializa B e D para a diagonal de A
 
        do ip=1,n
           b(ip)=a(ip,ip)
           d(ip)=b(ip)

!     Este vetor ira acumular termos da forma tapq como na equacao da Pag 458
!     do Numerical Recipes - eq (11.1.14) 
 
           z(ip)=0
        enddo
        nrot=0   
        do i=1,50     
           sm=0.

!     Soma dos elementos fora da diagonal

           do ip=1,n-1 
              do   iq=ip+1,n 
                   sm=sm+abs(a(ip,iq)) 
              enddo  
           enddo
           if(sm.eq.0.) goto 102
             if(i.lt.4)then
             tresh=0.2*sm/n**2
             else
             tresh=0
             endif
           do   ip=1,n-1
                do  iq=ip+1,n
                    g=100.*abs(a(ip,iq))

!      Apos quatro varreduras,pular a rotacao 
!      se o elemento fora da diagonal for pequeno.

            if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip)))            &
	 &      .and.(abs(d(iq))+g.eq.abs(d(iq)))) then                  
                a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh) then                     
				h=d(iq)-d(ip)
				if(abs(h)+g.eq.abs(h))then
				    t=a(ip,iq)/h   
                else

!      Equacao (11.1.10) da Pag-457 do Numerical Recipes

                   theta=0.5*h/a(ip,iq)
                   t=1./(abs(theta)+sqrt(1.+theta**2))
                   if(theta.lt.0.)t=-t
                   endif 
                   c=1./sqrt(1+t**2)
                   s=t*c
                   tau=s/(1.+c)
                   h=t*a(ip,iq)
                   z(ip)=z(ip)-h
                   z(iq)=z(iq)+h
                   d(ip)=d(ip)-h
                   d(iq)=d(iq)+h
                   a(ip,iq)=0

!       Caso de rotacoes 1<=j< p

                   do   j=1,ip-1
                   g=a(j,ip)
                   h=a(j,iq)
                   a(j,ip)=g-s*(h+g*tau)        
                   a(j,iq)=h+s*(g-h*tau)
                   enddo

!       Caso de rotacoes p < j < q

                   do   j=ip+1,iq-1    
                        g=a(ip,j)
                        h=a(j,iq)
                        a(ip,j)=g-s*(h+g*tau)
                        a(j,iq)=h+s*(g-h*tau)
                   enddo

!       Caso de rotacoes q < j <= n

                   do   j=iq+1,n
                        g=a(ip,j)
                        h=a(iq,j)
                        a(ip,j)=g-s*(h+g*tau) 
                        a(iq,j)=h+s*(g-h*tau)
                   enddo 
                   do   j=1,n
                        g=v(j,ip)
                        h=v(j,iq)
                        v(j,ip)=g-s*(h+g*tau)
                        v(j,iq)=h+s*(g-h*tau)
                   enddo    
                        nrot=nrot+1
           endif    
           enddo
           enddo
           do  ip=1,n
               b(ip)=b(ip)+z(ip)

!       Atualizar d com a soma de tapq e reinicializar z.
 
               d(ip)=b(ip)
               z(ip)=0
           ENDDO
       ENDDO 
!       pause 'muitas interacoes no JACOBI'

! A subrotina abaixo ordena os autovalores encontrados na Jacobi. 

102    CONTINUE
       Call EIGSRT(d,v,N)
       
!       Constroi matriz diagonalizada SD a partir de D
    do i=1,N
        do j=1,N
            SD(i,j)=0
        enddo
        SD(i,i)=D(i)
    enddo

       Return
       END SUBROUTINE jacobi 

    subroutine eigsrt(d,v,n) 
! Subrotina para ordenar os auovalores de Jacobi.(Ordem decrescente)
       
	   integer n
	   integer i,j,k
	   real * 8 d(n),v(n,n)
	   real * 8 p
	   do i=1,n-1
		 k=i
		 p=d(i)
		 do j=i+1,n
		   if (d(j).le.p) then
			   k=j        
			   p=d(j)
		   endif
		 enddo
		 if (k.ne.i) then
			 d(k)=d(i)
			 d(i)=p
			 do j=1,n
			   p=v(j,i)
			   v(j,i)=v(j,k)
			   v(j,k)=p
			 enddo
		 endif
	   enddo
	   return
	end subroutine eigsrt
    
    subroutine multMATRIX(MX1,I1,J1,MX2,I2,J2,MXF)
      
    integer i,j,k
    integer, intent(in) :: I1,I2,J1,J2  
    real(8), intent(in) :: MX1(I1,J1), MX2(I2,J2)
    real(8), intent(out) :: MXF(I1,J2)

!   Verifying that the matrices can be multiplied
    if (J1 .NE. I2) then
        print *, "Your matrices can not be multiplied, they do not have proper dimensions"
        print *, "Press enter to continue"
        read (*,*)
        stop
    endif

    MXF = 0.0
    
!   Multiplying the matrix M1(I1,J1) by the matrix M2(I2,J2)
    
    do i=1,I1
        do j=1,J2
            do k=1,I2
                MXF(i,j) = MXF(i,j) + (MX1(i,k)*MX2(k,j))
            enddo
        enddo
    enddo
    
    end subroutine multMATRIX
    
    subroutine simpson(XY,NPONTOS,H,S)
    
! Subrotina Simpson: Integra a funcao usando o metodo de Simpson

      integer npontos,i
      real*8 xy(1000),h,s
      
      s=0
      do 101 i=1,npontos-2,2
       s=s+xy(i)+4.*xy(i+1)+xy(i+2)
101   continue
      s=s*h/3.
      return
    end
    
end module