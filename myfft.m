function y = myfft(x)
L = length(x);
x1 = fft(x);
P2=abs(x1/L);
X1=P2;
X1(L/2:end)=P2(1:L/2+1);
X1(1:L/2)=P2(L/2+1:end);
% f = fs*(-L/2:L/2-1)/L;
% plotting power spectral:
% figure()
% plot(f,(abs(X1)).^2);
% xlabel('f(Hz)');
% ylabel('S(f)');
% title('problem 1 PART 4');
y = X1;
end