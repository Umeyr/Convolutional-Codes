u = 1;

actual = zeros(10,1);
actualuncoded = zeros(10,1);
theory = zeros(10,1);
pvalue = zeros(10,1);
snrvalue = zeros(10,1);

theorySNR = zeros(10,1);


for p = -0.3:-.1:-3
   
 SNR = 10;

  pvalue(u) = p;  
  snrvalue(u) = SNR;
 %Encoder part

h = 1000; %length of the input

m = randi([0 1],h,1);
m(h) = 0;
m(h-1) = 0;
c = zeros(h,3);

%uncoded transmisson
nuncoded = randn(h,1);

Esuncoded = sum((m(:,1).^2)) ;
Esuncoded = Esuncoded/(h);
 
nuncoded = nuncoded.*sqrt(Esuncoded/SNR);

yuncoded = m + nuncoded;

for i = 1:h
        if yuncoded(i,1) >= .5
            yuncoded(i,1 ) = 1;
        else
            yuncoded(i,1 ) = 0;
        end

end

actualuncoded(u) = sum(abs(m-yuncoded));


%%%%%%%%%%%%%%%%%%%%

state = zeros(1,3);

for i = 1:length(m)
   
    state(3) = state(2);
    state(2) = state(1);
    state(1) = m(i);
    
    c(i,1) = state(3);
    c(i,2) = state(1) + state(3);
    c(i,3) = sum(state);
end

c = rem(c,2);

p = 10^p;

cbsc = bsc(c,p); % Binary symmetric channel( this part used in question 2)
%cbsc = c;

%cbsc = 2*cbsc - 1;

%AWGN channel addition ( this part used in question 3)

n = randn(h,3);

Es = sum((cbsc(:,1).^2)) + sum((cbsc(:,2).^2)) + sum((cbsc(:,3).^2));
Es = Es/(3*h);

  
n = n.*sqrt(Es/SNR);

y = cbsc + n;

%y = awgn(cbsc,3*SNR,'measured');

%y = awgn(cbsc,10,'measured');
%%%%%%%%%%ENCODER is completed%%%%%%%%%%%%%%%%%%%%%%%%%%

%DECODER PART

crcvd = zeros(h,3);

for i = 1:h
    if cbsc(i,1) >= .5
        crcvd(i,1 ) = 1;
    end

    if cbsc(i,2) >= .5
        crcvd(i,2 ) = 1;
    end
    if cbsc(i,3) >= .5
        crcvd(i,3 ) = 1;
    end
end

%for BPSK in SDD only
%crcvd = 2*crcvd - 1;

%quantization
%two bit
% for i = 1:h
%     
%     if y(i,1) >= .5
%         y(i,1 ) = .75;
%     elseif y(i,1) >= 0
%         y(i,1 ) = .25;
%      elseif y(i,1) >= -0.5
%         y(i,1 ) = -0.25;
%     else
%         y(i,1 ) = -0.75;
%     end
% 
%     if y(i,2) >= .5
%         y(i,2 ) = .75;
%     elseif y(i,2) >= 0
%         y(i,2 ) = .25;
%      elseif y(i,2) >= -0.5
%         y(i,2 ) = -0.25;
%     else
%         y(i,2 ) = -0.75;
%     end
%     
%      if y(i,3) >= .5
%         y(i,3 ) = .75;
%     elseif y(i,3) >= 0
%         y(i,3 ) = .25;
%      elseif y(i,3) >= -0.5
%         y(i,3 ) = -0.25;
%     else
%         y(i,3 ) = -0.75;
%     end
%     
% end
%traceback in the trellis diagram

%three bit
% for i = 1:h
%     
%     if y(i,1) >= .75
%         y(i,1 ) = .875;
%     elseif y(i,1) >= 0.5
%         y(i,1 ) = 0.625;
%      elseif y(i,1) >= 0.25
%         y(i,1 ) = 0.375;
%      elseif y(i,1) >= 0
%         y(i,1 ) = 0.125;
%      elseif y(i,1) >= -0.25
%         y(i,1 ) = -0.125;
%      elseif y(i,1) >= -0.5
%         y(i,1 ) = -0.375;
%      elseif y(i,1) >= -0.75
%         y(i,1 ) = -0.625;
%     else
%         y(i,1 ) = -0.875;
%     end
%     
%     
%     if y(i,2) >= .75
%         y(i,2 ) = .875;
%     elseif y(i,2) >= 0.5
%         y(i,2 ) = 0.625;
%      elseif y(i,2) >= 0.25
%         y(i,2 ) = 0.375;
%      elseif y(i,2) >= 0
%         y(i,2 ) = 0.125;
%      elseif y(i,2) >= -0.25
%         y(i,2 ) = -0.125;
%      elseif y(i,2) >= -0.5
%         y(i,2 ) = -0.375;
%      elseif y(i,2) >= -0.75
%         y(i,2 ) = -0.625;
%     else
%         y(i,2 ) = -0.875;
%     end
%     
%     
%     if y(i,3) >= .75
%         y(i,3 ) = .875;
%     elseif y(i,3) >= 0.5
%         y(i,3 ) = 0.625;
%      elseif y(i,3) >= 0.25
%         y(i,3 ) = 0.375;
%      elseif y(i,3) >= 0
%         y(i,3 ) = 0.125;
%      elseif y(i,3) >= -0.25
%         y(i,3 ) = -0.125;
%      elseif y(i,3) >= -0.5
%         y(i,3 ) = -0.375;
%      elseif y(i,3) >= -0.75
%         y(i,3 ) = -0.625;
%     else
%         y(i,3 ) = -0.875;
%     end
% 
%     
% end



mintrace = zeros(1,1);
laststate = 0;

for t = 1:h/5
    
    trace = zeros( 6, 32);
    trace(1,1) = laststate;
    distance = zeros(1,32);
    distance = distance - 1;
    curstate = zeros(1,3);

    %crcvdmo = crcvd(5*t-4:5*t,:);
    crcvdmo = cbsc(5*t-4:5*t,:);
    
for k = 1:5
   
    for i = 1:2^(k-1)
        
    if trace(k, i) == 0
        
        %copy 
      trace(1:k, 2^(k-1)+i) = trace(1:k,i);
      distance(2^(k-1)+i) = distance(i);
      
      trace(k+1, i) = 0;
      
      %distance
      if distance(i) == -1
          distance(i) = 0;
      end
      curstate = [0 0 0];
      %distance(i) = distance(i) + sum(abs(crcvdmo(k,1:3)-curstate));
      distance(i) = distance(i) + sum((crcvdmo(k,1:3).*curstate));
      
      
      
      trace(k+1, 2^(k-1)+i) = 2;
      
      %distance
      if distance(i) == -1
          distance(i) = 0;
      end
      curstate = [1 0 0];
      %distance(2^(k-1)+i) = distance(2^(k-1)+i) + sum(abs(crcvdmo(k,1:3)-curstate));
     distance(2^(k-1)+i) = distance(2^(k-1)+i) + sum((crcvdmo(k,1:3).*curstate));
      
    elseif trace(k, i)== 1
        
         %copy 
         trace(1:k, 2^(k-1)+i) = trace(1:k,i);
         distance(2^(k-1)+i) = distance(i);
         
         trace(k+1, i) = 0;
         
         %distance
         if distance(i) == -1
          distance(i) = 0;
         end
        curstate = [0 1 1];
        %distance(i) = distance(i) + sum(abs(crcvdmo(k,1:3)-curstate));
        distance(i) = distance(i) + sum((crcvdmo(k,1:3).*curstate));
         
      
         trace(k+1, 2^(k-1)+i) = 2;
         
         %distance
         if distance(i) == -1
          distance(i) = 0;
          end
        curstate = [1 1 0];
        distance(2^(k-1)+i) = distance(2^(k-1)+i) + sum((crcvdmo(k,1:3).*curstate));
    
        
    elseif trace(k, i) == 2
        
        %copy 
          trace(1:k, 2^(k-1)+i) = trace(1:k,i);
          distance(2^(k-1)+i) = distance(i);
          
         trace(k+1, i) = 1;
         
         %distance
         if distance(i) == -1
          distance(i) = 0;
          end
         curstate = [1 1 0];
         distance(i) = distance(i) +  sum((crcvdmo(k,1:3).*curstate));
        
          
      
         trace(k+1, 2^(k-1)+i)  = 3;
         
         %distance
         if distance(i) == -1
          distance(i) = 0;
          end
         curstate = [0 1 1];
         distance(2^(k-1)+i) = distance(2^(k-1)+i) + sum((crcvdmo(k,1:3).*curstate));
    
     
    elseif trace(k, i) == 3
        
         %copy 
        trace(1:k, 2^(k-1)+i) = trace(1:k,i);
        distance(2^(k-1)+i) = distance(i);
        
         trace(k+1, i)= 1;
         
         %distance
         if distance(i) == -1
          distance(i) = 0;
         end
        curstate = [1 0 1];
        distance(i) = distance(i) + sum(abs(crcvdmo(k,1:3).*curstate));
        
         
      
        trace(k+1, 2^(k-1)+i)  = 3;
        
        %distance
        if distance(i) == -1
          distance(i) = 0;
        end
        curstate = [0 0 0];
        distance(2^(k-1)+i) = distance(2^(k-1)+i) + sum((crcvdmo(k,1:3).*curstate));
    
    end
           
    end
end

%[M,I] = min(distance);
[M,I] = max(distance);

mintrace = [mintrace;trace(2:6,I)];

laststate = trace(6,I);

end


received = zeros(h,1);

for i = 1:h
    
    if mintrace(i) == 0
        if mintrace(i+1) == 0
            received(i) = 0;
        elseif mintrace(i+1) == 2
            received(i) = 1;
        end
        
    elseif mintrace(i) == 1
        if mintrace(i+1) == 0
            received(i) = 0;
        elseif mintrace(i+1) == 2
            received(i) = 1;
        end
        
    elseif mintrace(i) == 2
        
        if mintrace(i+1) == 1
            received(i) = 0;
        elseif mintrace(i+1) == 3
            received(i) = 1;
        end
        
    elseif mintrace(i) == 3
        if mintrace(i+1) == 1
            received(i) = 0;
        elseif mintrace(i+1) == 3
            received(i) = 1;
        end
    end
end  


actual(u) = sum(abs(m-received));
theory(u) = 64*(p^3)*((1-p)^3)/(1-8*p*(1-p))^2;
theorySNR(u) = (exp(-2*SNR))/((1-2*exp(-2*SNR/3))^2);

u = u+1;
end

actual = actual/1000;

actual = log10(actual);

actualuncoded = actualuncoded/1000;

actualuncoded = log10(actualuncoded);

theory = log10(theory);

theorySNR = log10(theorySNR);



% figure(1);
% plot(pvalue, theory);
% title('Upper Bound for Bit Error Probability');
% xlabel('cross over probability p in dB');
% ylabel('Bit Error Probability in dB');

figure(1);
plot(pvalue, actual);
title('Experimental Bit Error Probability in SDD ');
xlabel('Eb/N0 in dB');
ylabel('Bit Error Probability in log scale');

% figure(2);hold on
% a1 = plot(pvalue, actual); M1 = 'SDD Result          ';
% a2 = plot(pvalue, actualuncoded); M2 = 'Uncoded Transmission';
% legend([a1,a2], [M1; M2]);
% title('Experimental Bit Error Probability');
% xlabel('Eb/N0 in dB');
% ylabel('Bit Error Probability in log scale');
% 
% figure(3);
% plot(pvalue, theorySNR);
% title('Upper Bound for Bit Error Probability for AWGN channel ');
% xlabel('Eb/N0 in dB');
% ylabel('Bit Error Probability in log scale');



