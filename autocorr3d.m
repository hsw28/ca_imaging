function f = autocorr3d(Data, minlag, maxlag, laginterval);

tic
% List of data arrays called 'Data'.

New_Data = [];

O = ones(size(Data,1), size(Data,2)); %np.ones((len(Data[0][0]),len(Data[0][1])))  Resolution of data arrays.
A = mean(Data,3);                        % Average 2D array
M = mean(A(:));                                 % Average of average 2D array
S = std(A);                                  % Standard deviation of average 2D array


New_Data = NaN(size(Data,1), size(Data,2), size(Data,3));
for i=1:(size(Data,3))
    New_Data(:,:,i) =  (Data(:,:,i)-(M.*O))./S; %New_Data.append((Data[i]-M*O)/S) # (X_t-mu)/sigma
    [x,y] = find(New_Data(:,:,i)==-inf);
    New_Data(x,y,i) = NaN;
end


%Autocorrelation for varying time lags.
R = [];
Count = minlag;
while Count<=maxlag+1 % Arbitrary choice for max lag time.
    Matrix_Multiply = NaN(size(Data,1), size(Data,2), size(Data,3));
    for j=1:(size(New_Data,3)-Count)
        Matrix_Multiply(:,:,j)=(New_Data(:,:,j).*New_Data(:,:,j+Count));
    end
    R(:,:,end+1) = [(nansum(Matrix_Multiply, 3))]./j;
    Count = Count+laginterval
end



Solution = [];
for k =1:(size(R,3))
  mean(R(:,:,k),'all');
    Solution(end+1)=mean(R(:,:,k),'all');
end

f = Solution;

toc

%t = [0.1*k for k in range(1,len(Solution)+1)]
figure
%plt.semilogy(t, Solution)
plot(Solution)
xlabel('Lag time')
ylabel('Matrix sum')
title('Field autocorrelation over time')
