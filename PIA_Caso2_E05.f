      program PIA_E05_C2
      !Este programa calcula la sol. num. de una e.d.p. (uxx+uyy=Cte.)
       !Con u=u(x,y)

      implicit none
      real*8:: l,h,dx,dy,alfa,b,b1,b2,c1,c2,d1,d2,e1,e2,s,la,p,h1,h2 !variables de paso
      real*8:: h3,h4,h5,h6,h7,h8
      real*8:: f, f1, f2, f3, f4, f5, f6, f7  !funciones analiticas
      integer*8:: i,j,ny,nx,k,k1   !contadores
      real*8, allocatable, dimension (:,:):: a, a1, a2, a3, a4 !matrices nxn
      real*8, allocatable, dimension (:,:):: x1, y1 !vectores nx1
      
      OPEN( 7, FILE= "MatrizNum.txt")
      OPEN( 8, FILE= "MatrizAna.txt")
      OPEN( 9, FILE= "Pia_Caso2.gpl")
      !scrip
      write(9,*)"set xlabel ""x"" font ""Times-Roman, 10"""
      write(9,*)"set ylabel ""y"" font ""Times-Roman, 10"""
      write(9,*)"set zlabel ""u"" font ""Times-Roman, 10"""
      write(9,*)"set title ""Pia caso 2"""
      write(9,*)"set xtics font ""Times-Roman, 10"""
      write(9,*)"set ytics font ""Times-Roman, 10"""
      write(9,*)"set ztics font ""Times-Roman, 10"""
      write(9,*)"splot ""MatrizNum.txt"", ""MatrizAna.txt"""

      write(*,*)"******************************************************"
      write(*,*)"Escriba los nodos en y (ny):"
      read(*,*) ny
      write(*,*)"Escriba los nodos en x (nx):"
      read(*,*) nx
      write(*,*)"Escriba la long. en x (l):"
      read(*,*) l
      write(*,*)"Escriba la long. en y (h):"
      read(*,*) h
      write(*,*)"Escriba el termino lamda:"
      read(*,*) la
      write(*,*)"Escriba las condiciones Dirichet:"
      write(*,*)"En u(0,y), u(l,y):"
      write(*,*)"Pd. estas se calcula automaticamento, porque,
     &  u(0,y)=u1(y) ^ u(l,y)=u2(y)"
      !read(*,*)c1, c2
      write(*,*)"En u(x,0), u(x,h):"
      read(*,*)d1, d2
      write(*,*)"Escriba el numero de iteraciones el metodo (k)
     & si k>>1 las soluciones se parecen mas"
      read(*,*) k1
      write(*,*)"******************************************************"

      !fijacion de parametros
      p=4*atan(1.0)
      dx= l/nx
      dy= h/ny
      b=  2*(dx**2 + dy**2)
      s= ((dx*dy)**2)*la
      b1= dx**2
      b2= dy**2
      
      allocate( a(ny+1,nx+1) )
      allocate( a1(ny+1,nx+1) )
      allocate( x1(nx+1,1) )
      allocate( y1(ny+1,1) )

      !Llenado de las condiciones de frontera
      Do i=1, ny+1
      Do j=1, nx+1
      
      c1= i*dy   !incremento con los que se evalua u(x,0)^ u(x,l)
      !c2= (i)*dy
      a(i,j)=0
      a(i,1)= -1*f4(c1,la,h)  !condiciones en x
      a(i,nx+1)= -1*f4(c1,la,h)
      
      a(1,j)= d1   !condiciones en y
      a(ny+1,j)=d2


      end do
      end do

      !Llenado de vectores
      Do i=1, ny+1
      y1(i,1)=i*dy
      End do

      Do i=1, nx+1
      x1(i,1)=i*dx
      End do

      !Matriz de los puntos interiores
      Do k=1,k1
      Do i=2,ny
      Do j=2,nx        !cambie el b1 y b2
                !formula de DFM
      a(i,j)=(b2/b)*(a(i+1,j)+a(i-1,j))+(b1/b)*(a(i,j+1)+a(i,j-1)) + s/b

      end do
      end do
      end do
      
      !Matriz de la sol. ana.
      Do i=1, ny+1
      Do j=1, nx+1

      e1=j*dx
      e2=i*dy
      a1(i,j)= f(e1,e2,la,h,l)
      
       !incremento con los que se evalua u(x,0)^ u(x,l)

      a1(i,1)= -1*f4(e2,la,h)  !condiciones en x
      a1(i,nx+1)= -1*f4(e2,la,h)

      a1(1,j)= d1   !condiciones en y
      a1(ny+1,j)= d2

      end do
      end do


      !Escritura de matriz y vectores
      write(*,*)"Matriz de resultados:"
      write(*,*)"Numericos:"
      Do i= 1,ny+1
      Do j=1,nx+1
      !write(*,*) x1(j,1), y1(i,1), a(i,j)
      write(7,*) x1(j,1), y1(i,1), a(i,j)
      end do
      end do
      
      write(*,*)"Matriz de resultados:"
      write(*,*)"Analiticos:"
      Do i= 1,ny+1
      Do j= 1,nx+1
      !write(*,*) x1(j,1), y1(i,1), a1(i,j)
      write(8,*) x1(j,1), y1(i,1), a1(i,j)
      end do
      end do

      close(7)
      close(8)
      close(9)
      call system("Pia_Caso2.gpl")

      write(*,*)"Aqui acaba el programa:)"
      pause
      stop
      continue
      end
       
       
       function f1(x,y,h,l)
       implicit none
       real*8:: f,f1,f2,f3,x,y,la,h,l,p,h1,h2,h3,h4
       !Parametros sol. ana.
       p=4*atan(1.0)
       h2=1/( sinh(p*l/h) )

       f1= ( sinh(p*x/h) + sinh(p*(l-x)) )*( sin(p*y/h) )*h2

       return
       end
       
       function f2(x,y,h,l)
       implicit none
       real*8::  f,f1,f2,f3,x,y,la,h,l,p,h1,h2,h3,h4
       !Parametros sol. ana.
       p=4*atan(1.0)
       h3=1/( 9*sinh(p*3*l/h)  )

       f2= ( sinh(p*3*x/h) + sinh(p*3*(l-x)) )*( sin(p*3*y/h) )*h3

       return
       end

       function f3(x,y,h,l)
       implicit none
       real*8::  f,f1,f2,f3,x,y,la,h,l,p,h1,h2,h3,h4
       !Parametros sol. ana.
       p=4*atan(1.0)
       h4=1/( 125*sinh(p*5*l/h) )

       f3= ( sinh(p*5*x/h) + sinh(p*5*(l-x)) )*( sin(p*5*y/h) )*h4

       return
       end
       
       function f5(x,y,h,l)
       implicit none
       real*8::  f,f1,f2,f3,f5,f6,x,y,la,h,l,p,h1,h2,h3,h4,h5,h6
       !Parametros sol. ana.
       p=4*atan(1.0)
       h6=1/(343*sinh(p*7*l/h) )

       f5= ( sinh(p*7*x/h) + sinh(p*7*(l-x)) )*( sin(p*7*y/h) )*h6

       return
       end
       
       function f6(x,y,h,l)
       implicit none
       real*8::  f,f1,f2,f3,f5,f6,x,y,la,h,l,p,h1,h2,h3,h4,h5,h6,h7
       !Parametros sol. ana.
       p=4*atan(1.0)
       h7=1/( 729*sinh(p*9*l/h) )

       f6= ( sinh(p*9*x/h) + sinh(p*9*(l-x)) )*( sin(p*9*y/h) )*h7

       return
       end
       
       function f7(x,y,h,l)
       implicit none
       real*8::f,f1,f2,f3,f5,f6,f7,x,y,la,h,l,p,h1,h2,h3,h4,h5,h6,h7,h8
       !Parametros sol. ana.
       p=4*atan(1.0)
       h8=1/( 1331*sinh(p*11*l/h) )

       f7= ( sinh(p*11*x/h) + sinh(p*11*(l-x)) )*( sin(p*11*y/h) )*h8

       return
       end
       
       !splot "MatrizAna.txt", "MatrizNum.txt"

       !condiciones de borde
       function f4(y,la,h)
       implicit none
       real*8:: f4,y,la,h,h1
       !Parametros sol. ana.
       h1=0.5*(la*h**2)

       f4=h1*( (y/h)-(y/h)**2 )

       return
       end
      
              !Sol. analitica
       function f(x,y,la,h,l)
       implicit none
       real*8:: f,f1,f2,f3,f4,f5,f6,f7,x,y,la,h,l,p,h1,h2,h3,h4,h5,h6
       p=4*atan(1.0)
       h5=(4*la*h**2)/((p**3))

       f= -1*f4(y,la,h)-(h5)*f1(x,y,h,l) -h5*f2(x,y,h,l)
     & -h5*f3(x,y,h,l) - h5*f5(x,y,h,l) - h5*f6(x,y,h,l)+h5*f7(x,y,h,l)
       return
       end
      
      
      
      


