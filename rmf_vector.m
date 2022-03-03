function [ret] = rmf_vector(kmu) 

    kap = sqrt(sum(kmu.^2));
    mu = kmu./kap;
    m = length(mu);
    if (kap == 0)
        u = randn(length(kmu),1);
        u = u./sqrt(sum(u.^2));
    end
    if (kap > 0)
        if (m == 1)
            u = (-1)^(binornd(1,1/(1 + exp(2 * kap * mu))));
        end
        if (m > 1) 
            W = rW(kap, m);
            V = randn(m - 1,1);
            V = V/sqrt(sum(V.^2));
            x = [(1 - W^2)^0.5*V' W];
            u = [null(mu) mu'] * x';
        end
    end
    ret = u;
end

function [W] = rW(kap, m)
  W = 1.0;
  b=( -2.0*(kap) + sqrt(  4*kap^2+(m-1.0)^2 ) )/(m-1.0) ;
  x0=(1.0-b)/(1.0+b) ;
  c= (kap)*x0 +(m-1.0)*log(1.0-(x0^2)) ;
  done=false;

  while( done==0)
    Z = betarnd( (m-1.0)/2.0, (m-1.0)/2.0 );
    W=( 1-(1+b)*Z)/(1.0-(1.0-b)*Z) ;
    U=rand;
    if( (kap)*(W)+(m-1)*log(1-x0*(W))-c  > log(U) )
        done=true;
    end
  end
end