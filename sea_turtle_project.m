%sea turtle population dynamics analyzed with Lefkovitch matrices
%initialize the population number vector (life stages 1-7)
n = [100;150;40;20;15;10;5];
%initialize the population matrix
L = [0,0,0,0,127,6,95;.68, .74,0,0,0,0,0;0,0.05,.69,0,0,0,0;0,0,.01,.67,0,0,0;...
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
        e(j,k) = L(j,k)*W(j,:)*V(:,j)/lambda;
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
figure
subplot(2,2,1)
bar(Ps)
title('Elasticity of lambda vs. Changes in Ps')
xlabel('Life Stage Number')
ylabel('Elasticity of lambda')
subplot(2,2,2)
bar(Gs)
title('Elasticity of lambda vs. Changes in Gs')
xlabel('Life Stage Number')
ylabel('Elasticity of lambda')
subplot(2,2,3)
bar(Fs)
title('Elasticity of Principal Eigenvalue vs. Changes Fs')
xlabel('Life Stage Number')
ylabel('Elasticity of lambda')
saveas(gcf, 'turtleBars.pdf')
for i = 1:100
    ni(i,:) = L^i*n;
    n1(i) = ni(i,1);
    n2(i) = ni(i,2);
    n3(i) = ni(i,3);
    n4(i) = ni(i,4);
    n5(i) = ni(i,5);
    n6(i) = ni(i,6);
    n7(i) = ni(i,7);
end
t = 1:100;
figure
plot(t,n1,t,n2,t,n3,t,n4,t,n5,t,n6,t,n7)
title('Number of Turtles in Each Life Stage From t = 0 to 1 = 100')
xlabel('time')
ylabel('number of turtles in each life stage')
legend('Hatchlings','Small Juveniles','Large Juveniles','Sub-Adults','Novice Breeders','1st-year Remigrants','Mature Breeders')

figure
semilogy(t,n1,t,n2,t,n3,t,n4,t,n5,t,n6,t,n7)
title('Number of Turtles in Each Life Stage From t = 0 to 1 = 100 (logarithmic plot)')
xlabel('time')
ylabel('number of turtles in each life stage')
legend('Hatchlings','Small Juveniles','Large Juveniles','Sub-Adults','Novice Breeders','1st-year Remigrants','Mature Breeders')

x = [900;310;443;220;875;123;501];    
for i = 1:100
    xi(i,:) = L^i*x;
    x1(i) = xi(i,1);
    x2(i) = xi(i,2);
    x3(i) = xi(i,3);
    x4(i) = xi(i,4);
    x5(i) = xi(i,5);
    x6(i) = xi(i,6);
    x7(i) = xi(i,7);
end
figure
plot(t,x1,t,x2,t,x3,t,x4,t,x5,t,x6,t,x7)
title('Number of Turtles in Each Life Stage From t = 0 to 1 = 100')
xlabel('time')
ylabel('number of turtles in each life stage')
legend('Hatchlings','Small Juveniles','Large Juveniles','Sub-Adults','Novice Breeders','1st-year Remigrants','Mature Breeders')        
        

        
        
        
        
