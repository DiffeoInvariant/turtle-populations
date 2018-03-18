%sea turtle population dynamics analyzed with Lefkovitch matrices

%initialize the population number vector (life stages 1-7)
n = [100;150;40;20;15;10;5];
%initialize the population matrix
L = [0,0,0,0,127,6,95;.68, .74,0,0,0,0,0;0,0.05,.67,0,0,0,0;0,0,.01,.67,0,0,0;...
    0,0,0,.05,0,0,0;0,0,0,0,.82,0,0;0,0,0,0,0,.79,.83];
%Compute the matrix of eigenvectors (V) and diagonal matrix of eigenvalues (D)
%and the matrix of left eigenvectors, W
[V,D,W] = eig(L);
%make the eigenvectors real-only
[a,b] = size(V);
V = real(V);
W = real(W);
%find the largest eigenvalue lambda, indexed by i
[lambda, i] = max(D(:));
[irow, icol] = ind2sub(size(D),i);
%find the eigenvector corresponding to it
vmax = V(:,icol);
%normalize
vmax = vmax/sum(vmax);
%do the same for the left eigenvector
wmax = W(irow,:);
wmax = wmax/sum(wmax);
%create the sensitivity matrix
s = zeros(7);
for j = 1:7
    for k = 1:7
        s(j,k) = wmax(j)*vmax(k)/(wmax*vmax);
    end
end
%create the elasticity matrix
e = zeros(7);
for j = 1:7
    for k = 1:7
        e(j,k) = L(j,k)*s(j,k)/lambda;
    end
end
%create vectors of Ps, Gs, and Fs sensitivities in each stage.
Ps = zeros(7,1);
Gs = zeros(6,1);
Fs = zeros(7,1);
for j = 1:7
        Ps(j) = s(j,j);
        if j ~= 7
            Gs(j) = s(j+1,j);
        end
        Fs(j) = s(1,j);
end
subplot(2,2,1)
bar(Ps)
subplot(2,2,2)
bar(Gs)
subplot(2,2,3)
bar(Fs)

        
        
        
        