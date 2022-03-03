function [ ret ] = rmf_matrix(M) 
% Generates a random orthonormal matrix from the von-Mises-Fisher
% distribution with mode M
    if (size(M,2) == 1)
        %1d case
        X = rmf_vector(M');
    else
        [UU,DD,VV] = svd(M);
        diagDD = diag(DD);
        H = UU * DD;
        m = size(H,1);
        R = size(H,2);
        cmet = false;
        trm = false;
        rej = 0;
        while (~cmet)
            U = zeros(m, R);
            U(:, 1) = rmf_vector(H(:,1)');
            lr = 0;
            for (j = linspace(2,R,R-1))
%                 N = null(U(:,linspace(1,j-1,j-1))');
                ZZ = U(:,linspace(1,j-1,j-1))';
                [Q,~,~] = qr(ZZ');
                N = Q(:,j:m);
                x = rmf_vector( (N'*H(:,j))' );
                U(:,j) = N * x;
                if (diagDD(j) > 0)
                  xn = sqrt(sum((N'*H(:,j)).^2));
                  xd = sqrt(sum(H(:,j).^2));
                  lbr = log(besseli( 0.5 * (m - j - 1), xn, 1)) - log( besseli(0.5 * (m - j - 1), xd, 1));
%                  if (is.na(lbr))
%                     lbr = 0.5 * (log(xd) - log(xn));
%                   end
                  lr = lr + lbr + (xn - xd) + 0.5 * (m - j - 1) * (log(xd) - log(xn));
               end
            end
            cmet = (log(rand(1)) < lr);
            rej = rej + (1 - 1 * cmet); %number of rejected samples
%             disp(rej)
            
            %if too many rejections, just return mode
            if (rej >= 1e2)
                cmet = true;
                trm = true;
            end
        end
        if (~trm)
            X = U * VV';
        else
            X = M;
        end
    end
    ret = X;
end