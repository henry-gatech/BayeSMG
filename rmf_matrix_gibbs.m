function [ ret ] = rmf_matrix_gibbs(M,X,rscol) 
% Generates a random orthonormal matrix from the von-Mises-Fisher
% distribution with mode M via Gibbs, initialized at X. rscol is the
% number of columns to be updated simultaneously

    if (rscol < 0)
        rscol = 2;
    end
    [sMu,sMd,sMv] = svd(M);
    H = sMu * sMd;
    Y = X * sMv;
    m = size(H,1);
    R = size(H,2);
    
%     succ_flg = false; %success flag
%     tol = 1e-1;
    
%     while (~succ_flg)
        for (iter = 1:round(R/rscol))
            r = randsample(R,rscol);
            idxmr = 1:size(Y,2);
            idxmr = idxmr(~ismember(idxmr,r));
    %         N = null(Y(:,idxmr)');
            ZZ = Y(:,idxmr)';
            [Q,~,~] = qr(ZZ');
            N = Q(:,(length(idxmr)+1):size(Y,1));
            y = rmf_matrix(N'*H(:,r));
            Y(:,r) = N * y;
        end
        ret = Y * sMv';

%         %check orthogonality
%         ck = ret'*ret;
%         if (norm(ck-eye(R)) < tol)
%             disp('Numerical error in Matrix Gibbs... returning mode.')
%             ret = M;
%         end
        
%     end
    
end