function [Ar,Br,Cr,Dr,HSVs] =   ERA(YY,m,n,numInputs,numOutputs,r)

for i=1:numOutputs
     for j=1:numInputs
         Dr(i,j) = YY(i,j,1);
         Y(i,j,:) = YY(i,j,2:end);
     end
end

assert(length(Y(:,1,1))==numOutputs);
assert(length(Y(1,:,1))==numInputs);
assert(length(Y(1,1,:))>=m+n);

for i=1:m
    for j=1:n
        for Q=1:numOutputs
            for P=1:numInputs
                H(numOutputs*i-numOutputs+Q,numInputs*j-numInputs+P) = Y(Q,P,i+j-1);
                H2(numOutputs*i-numOutputs+Q,numInputs*j-numInputs+P) = Y(Q,P,i+j);
            end
        end
    end
end

[U,S,V] = svd(H,'econ');
Sigma = S(1:r,1:r);
Ur = U(:,1:r);
Vr = V(:,1:r);
Ar = Sigma^(-.5)*Ur'*H2*Vr*Sigma^(-.5);
Br = Sigma^(-.5)*Ur'*H(:,1:numInputs);
Cr = H(1:numOutputs,:)*Vr*Sigma^(-.5);
HSVs = diag(S);