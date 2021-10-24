%f = @(x) 100 * (x {2} - x {1} .^ 2) .^ 2 + (1 - x {1}) .^ 2

function [xmin, fmin, maxit] = simplex(f, A, alpha, beta, gamma, maxit)
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
        fmin = min(feval);
        fg = max(feval(feval ~= max(feval)));
        %feval
        k = 0;
        do 
            k++;
        until (feval(k) == fmax)
        m = 0;
        do 
            m++;
        until(feval(m) == fmin)
        n = 0;
        do 
            n++;
        until(feval(n) == fg)
        xh = A{k};
        xl = A{m};
        xg = A(n);
        xc = centroide(A, k);
        xr = reflejado(A, xc, k);
        A{k} = xr;
        if(f(xr) < f(A{m}))
            A{k} = expansion(gamma, xc, A, k);
        elseif(f(xr) >= f(A{k}))
            A{k} = contraccion1(beta, xc, A, k);
        elseif(f(A{n}) < f(xr) && f(xr) < f(A{k}))
            A{k} = contraccion2(beta, xc, A, k);
        endif
        contador += 1;
    until (contador > maxit)
    A{k}
    f(A{k})
    maxit
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
            suma += 1/(N - 1) * x{i}{coord};
        endif
    endfor
endfunction

function [xr] = reflejado(x, xc, k)
    N = ndims(x);
    for j = 1:N
        xr{j} = 2 * xc{1}{j} - x{k}{j};
    endfor
endfunction

function [xnew] = expansion(gamma, xc, x, k)
    N = ndims(x);
    for j = 1:N
        xnew{j} = (1 + gamma) * xc{1}{j} - gamma * x{k}{j};
    endfor
endfunction

function [xnew] = contraccion1(beta, xc, x, k)
    N = ndims(x);
    for j = 1:N
        xnew{j} = (1 - beta) * xc{1}{j} + beta * x{k}{j};
    endfor
endfunction

function [xnew] = contraccion2(beta, xc, x, k)
    N = ndims(xc);
    for j = 1:N
        xnew{j} = (1 + beta) * xc{1}{j} - beta * x{k}{j};
    endfor
endfunction
            
            
