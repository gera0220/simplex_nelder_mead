%f = @(x) 100 * (x {2} - x {1} .^ 2) .^ 2 + (1 - x {1}) .^ 2

function [xmin, fmin] = simplex(f, A, alpha, maxit)
    feval = [];
    contador = 1;
    N = ndims(A);
    delta = alpha * [((sqrt(N + 1) + N - 1) / (N * sqrt(2))), ((sqrt(N + 1) - 1) / (N * sqrt(2)))];
    for i = 2:(N + 1)
        for j = 1:N
            if (i - 1 == j)
                A{i}{j} = A{1}{j} + delta(1);
            else
                A{i}{j} = A{1}{j} + delta(2);
            endif
        endfor
    endfor
    do
        %printf("Este es el %d simplex: ", contador), disp(A)
        for i = 1:length(A)
            feval(i) = f(A{i});
        endfor
        fmax = max(feval);
        if(exist("k") && max(feval) == f(A{k}))
            fmax = max(feval(feval ~= max(feval)));
        endif
        k = 0;
        do 
            k++;
        until (feval(k) == fmax)
        xc = centroide(A, k);
        A{k} = nuevo_punto(A, xc, k);
        contador += 1;
    until (contador > maxit)
    fmin = min(feval);
    k = 0;
    do
        k++;
    until(feval(k) == fmin)
    xmin = A{k};
    fmin = min(feval);
endfunction

function [xc] = centroide(x, k)
    N = ndims(x);
    for j = 1:N
        xc{1}{j} = suma_excepcion(x, j, k);
    endfor
endfunction

function [suma] = suma_excepcion(x, coord, k)
    suma = 0;
    N = length(x);
    for i = 1 : N
        if(k != i)
            suma += 0.5 * x{i}{coord};
        endif
    endfor
endfunction

function [xnew] = nuevo_punto(x, xc, k)
    N = ndims(x);
    for j = 1:N
        xnew{j} = 2 * xc{1}{j} - x{k}{j};
    endfor
endfunction


            
            
