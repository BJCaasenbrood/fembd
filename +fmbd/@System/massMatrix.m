function Mq = massMatrix(s, q)
% Mass matrix M(q) (constant for rigid bodies).
    Mq = zeros(s.nv);
    for i = 1:numel(s.units)
        B = s.units(i);  c = B.cols;
        if B.elastic
            Mq(c, c) = elasticTerms(B, q.R(:,:,i), q.e(B.ecols));
        else
            Mq(c, c) = blkdiag(B.m*eye(3), B.J);
        end
    end
end
