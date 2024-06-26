clear all; close all;
 
load('N50M100_30000');
X_50 = X; U_50 = U; U_ref_50 = U_ref;

figure;
subplot(1,2,1); hold on
scatter(X_50,U_50(end,:),'r'); scatter(X_50,U_ref_50,'b');
subplot(1,2,2);
scatter(X_50,abs(U_50(end,:)-U_ref_50));
title = ('N=50');

load('N100M100_30000');
X_100 = X; U_100 = U; U_ref_100 = U_ref;

figure;
subplot(1,2,1); hold on
scatter(X_100,U_100(end,:),'r'); scatter(X_100,U_ref_100,'b');
subplot(1,2,2);
scatter(X_100,abs(U_100(end,:)-U_ref_100));
title = ('N=100');