!defines real and integer kinds
module kinddef
    integer, parameter :: real_kind = 8
    !integer, parameter :: integer_kind = 4
endmodule


module adaptive_gauss_parameter
    use kinddef
    real(real_kind), dimension(*), parameter :: c = [  0.00600374098975726, 0.03136330379964697, 0.07589670829478634,&
                                                                            0.13779113431991491, 0.21451391369573058, 0.30292432646121836,&
                                                                            0.39940295300128270, 0.49999999999999994, 0.60059704699871719,&
                                                                            0.69707567353878153, 0.78548608630426942, 0.86220886568008526,&
                                                                            0.92410329170521366, 0.96863669620035286, 0.99399625901024291]
    real(real_kind), dimension(*,*), parameter ::  b = transpose(reshape([&
	                                                                            [ 0.015376620998058, 0.021474027715508, 0.000000000000000 ],&
	                                                                            [ 0.035183023744054, 0.014373156671397, 0.107157609486217 ],&
	                                                                            [ 0.053579610233586, 0.092599215594225, 0.000000000000000 ],&
	                                                                            [ 0.069785338963077, 0.011827744638032, 0.031130901929808 ],&
	                                                                            [ 0.083134602908497, 0.158470033694399, 0.000000000000000 ],&
	                                                                            [ 0.093080500007781, 0.003842920901955, 0.361711488583977 ],&
	                                                                            [ 0.099215742663556, 0.197412900561057, 0.000000000000000 ],&
	                                                                            [ 0.101289120962780, 0.000000000000000, 0.000000000000000 ],&
	                                                                            [ 0.099215742663556, 0.197412902276762, 0.000000000000000 ],&
	                                                                            [ 0.093080500007781, 0.003842917758458, 0.361711488583963 ],&
	                                                                            [ 0.083134602908497, 0.158470037700532, 0.000000000000000 ],&
	                                                                            [ 0.069785338963077, 0.011827740544188, 0.031130901929822 ],&
	                                                                            [ 0.053579610233586, 0.092599218962002, 0.000000000000000 ],&
	                                                                            [ 0.035183023744054, 0.014373154619738, 0.107157609486213 ],&
	                                                                            [ 0.015376620998058, 0.021474028361748, 0.000000000000000 ]&
                                                                            ], [3,15]))
endmodule


module Quadratur
        use kinddef
    interface
        elemental function func(x) result(res)
            use kinddef
            real(real_kind), intent(in) :: x
            real(real_kind) :: res
        endfunction
    endinterface
    
    contains
        !integrate one-dimensional function F from alpha to beta in N Segments with points c and weights b
        pure function quadrature_integrate(f, alpha, beta, n, b, c) result(intapprox)
            real(real_kind), intent(in) :: alpha, beta              !start and end of the interval
            real(real_kind) :: intapprox                                !result
            real(real_kind), dimension(:), intent(in) :: b, c    !weights and points
            integer, intent(in) :: n                                       !number of Segments in the decomposition of the interval
            procedure(func) :: f                                           !the function to integrate
            real(real_kind) :: h
            h = (beta - alpha) / N
            intapprox = sum(f( (/((i+c)*h+alpha, i=0,N-1)/) )* reshape(spread(b,2,N),(/30/)) * h)
        endfunction
        
        !integrate one-dimensional function F from alpha to beta with an adaptive quadrature with a relative error tolerance of tol
        pure function numint(f, alpha, beta, tol) result(approx)
        use adaptive_gauss_parameter, only : b, c   !import weights and absicasses of the 15-point Gauss-Legendre quadrature and weights for a 14 and 7 point quadrature with the same absicasses for error approximation
            real(real_kind), intent(in) :: alpha, beta                                                              !start and end of the interval
            real(real_kind), dimension(:), allocatable :: abs_approx , error, h, left, diff, diff_2
            real(real_kind), dimension(:,:), allocatable :: values                                             !evaluate the function once for all three quadrature rules and save it here
            real(real_kind) :: approx                                                                                     !result
            real(real_kind), intent(in) :: tol                                                                            !relative error tolerance
            integer :: maxerror
            procedure(func) :: f
            h = [(beta - alpha)]
            left = [alpha]
666      values = f( spread(c,2,size(h))*transpose(spread(h,2,size(c)))+transpose(spread(left,2,size(c))) )
            approx = sum(values * spread(b(:,1),2,size(h)) * transpose(spread(h,2,size(c))))
            abs_approx = sum(abs(values)* spread(b(:,1),2,size(h)) * transpose(spread(h,2,size(c))),1)
            diff = sum(values* spread(b(:,1) - b(:,2),2,size(h)) * transpose(spread(h,2,size(c))),1)
            diff_2 = sum(values* spread(b(:,1) - b(:,3),2,size(h)) * transpose(spread(h,2,size(c))),1)
            error = (diff / diff_2)**2
            if(sum(error) .ge. tol * sum(abs_approx)) then
                maxerror = maxloc(error, 1)
                h = [h(:maxerror - 1), h(maxerror) / 2._8, h(maxerror) / 2._8, h(maxerror + 1:)]
                left = [left(:maxerror), left(maxerror) + h(maxerror), left(maxerror + 1:)]
                goto 666
            endif
        endfunction
        
endmodule