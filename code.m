clc;
clear all;
close all;
n=importdata('ECG.txt');
x=n';
        
N = length (x);  % Signal length
fs = N/10;  % Sampling rate 
t = [0:N-1]/fs;  % time index

figure(1)
plot(x)
xlabel('Samples');ylabel('Amplitude');title('Input ECG Signal')
xlim([0 length(x)])

findpeaks(x)
xlabel('Samples')
ylabel('Amplitude')
title('Detecting Saturated Peaks')

x = x - mean (x );    % cancel DC conponents
x = x/ max( abs(x )); % normalize to one

figure(2)
plot(t,x)
xlabel('second');ylabel('Volts');title(' ECG Signal after cancellation DC drift and normalization')

% LPF (1-z^-6)^2/(1-z^-1)^2
b=[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a=[1 -2 1];


h_LP=filter(b,a,[1 zeros(1,12)]); % transfer function of LPF

x2 = conv (x ,h_LP);
%x2 = x2 (6+[1: N]); %cancle delay
x2 = x2/ max( abs(x2 )); % normalize , for convenience .

figure(3)
plot([0:length(x2)-1]/fs,x2)
xlabel('second');ylabel('Volts');title(' ECG Signal after LPF')
xlim([0 max(t)])
 
 
% HPF = Allpass-(Lowpass) = z^-16-[(1-z^-32)/(1-z^-1)]
b = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
a = [1 -1];

h_HP=filter(b,a,[1 zeros(1,32)]); % impulse response iof HPF

x3 = conv (x2 ,h_HP);
%x3 = x3 (16+[1: N]); %cancle delay
x3 = x3/ max( abs(x3 ));

figure(4)
plot([0:length(x3)-1]/fs,x3)
xlabel('second');ylabel('Volts');title(' ECG Signal after HPF')
xlim([0 max(t)])

xx=max(x);
possreg =(x > 0.60*xx)';

leftp = find(diff([0 possreg'])==1);
rightp = find(diff([possreg' 0])==-1);

for i=1:length(leftp)
    [Rvalue(i) Rloc(i)] = max( x(leftp(i):rightp(i)) );
    Rloc(i) = Rloc(i)-1+leftp(i); % add offset
end

Rloc=Rloc(find(Rloc~=0));
p=length(Rloc);
Rloc1=Rloc(1:(p));
p1=length(Rvalue)
Rvalue1=Rvalue(1:(p1));

figure(5)
plot ([0:length(x2)-1]/fs,x2,[0:length(x)-1]/fs,x,'r' )
legend('original signal', 'filtered signal');title('Original signal superimposed over filtered signal(notice the shift)')
xlabel('second');ylabel('Volts');
xlim([0,10]);

figure(6)
plot ([0:length(x)-1]/fs,x, Rloc1/fs ,Rvalue1 , 'r^');
legend('ECG','R peaks');title('ECG peak detection')
xlabel('second');ylabel('Volts');
xlim([0,10]);



xy=max(x);
possregg =(x > 0.60*xy)';

leftpp = find(diff([0 possregg'])==1);
rightpp = find(diff([possregg' 0])==-1);

for i=1:length(leftpp)
    [Rvaluep(i) Rlocp(i)] = max( x(leftpp(i):rightpp(i)) );
    Rlocp(i) = Rlocp(i)-1+leftpp(i); % add offset
end

Rlocp=Rlocp(find(Rlocp~=0));
pp=length(Rlocp);
Rloc1p=Rlocp(1:(pp));
p1p=length(Rvaluep)
Rvalue1p=Rvaluep(1:(p1p));


select=0;

if length(Rloc1p)==length(Rloc1)
    
select=select+1;    

figure(7)
plot (t,x, t(Rloc1p) ,Rvalue1p , 'r^');
legend('ECG','R peaks');title('ECG peak detection')
xlabel('second');ylabel('Volts');
xlim([0,10]);

else 
    
% g=Rloc1p-Rloc1;
% tdelay= t(Rloc1)-t(Rloc1p)
% del=1:length(tdelay)
% count=0
% for i=1:length(tdelay)
%    if tdelay(i)<0.0250   %%thsi section is for finding delay
%     del(i)=tdelay(i)
%     count=count+1;
%     else
%         del(i)=0
%     end
% end
% 
% avgdelay=sum(del)/count


figure(7)
plot ([0:length(x2)-1]/fs,x2, Rloc1/fs,Rvalue1 , 'r^');
legend('ECG','R peaks');title('ECG peak detection')
xlabel('second');ylabel('Volts');
xlim([0,10]);

end
 RLOCATION=[]; 
 RVALUE=[]
 xecg=[];
if select==1
RLOCATION=Rloc1p;
RVALUE=Rvalue1p;
xecg=x;
else
RLOCATION=Rloc1;
RVALUE=Rvalue1;
xecg=x2;
end
    
    

np=importdata('PPG1.txt');
xp=np';
        
Np = length (xp);  % Signal length
fsp = Np/10;  % Sampling rate 
tp = [0:Np-1]/fs;  % time index

figure(8)
plot(xp)
xlabel('Samples');ylabel('Amplitude');title('Input PPG Signal')
xlim([0 length(xp)])

xp = xp - mean (xp );    % cancel DC conponents
xp = xp/ max( abs(xp )); % normalize to one

figure(9)
plot(tp,xp)
xlabel('second');ylabel('Volts');title(' PPG Signal after cancellation DC drift and normalization')

% LPF (1-z^-6)^2/(1-z^-1)^2
b=[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a=[1 -2 1];


h_LP=filter(b,a,[1 zeros(1,12)]); % transfer function of LPF

x2p = conv (xp ,h_LP);
%x2 = x2 (6+[1: N]); %cancle delay
x2p = x2p/ max( abs(x2p )); % normalize , for convenience .

figure(10)
plot([0:length(x2p)-1]/fs,x2p)
xlabel('second');ylabel('Volts');title(' PPG Signal after LPF')
xlim([0 max(t)])
 
figure(11)
plot([0:length(x2p)-1]/fs,x2p, tp,xp,'r')
xlabel('second');ylabel('Volts');title(' PPG Signal Superimposition')
legend('Filtered PPG signal','Raw signal')
xlim([0 max(t)])
%  
% sx=diff(x2p)
% plot([0:length(x2p)-1]/fs,x2p, [0:length(sx)-1]/fs,sx,'r')
  
x_loc=[];
y_val=[];
x_max_loc=[];
y_max_val=[];
S_PEAK_LOC=[];
S_PEAK_VAL=[];
T_PEAK_LOC=[];
T_PEAK_VAL=[];
P_PEAK_LOC=[];
P_PEAK_VAL=[];
Q_PEAK_LOC=[];
Q_PEAK_VAL=[];
ht=(RLOCATION(2)-RLOCATION(1))/4;
for k =1:(length(RLOCATION)-1)
x_m=x2p(RLOCATION(k):RLOCATION(k+1));
x_e=xecg(RLOCATION(k):RLOCATION(k)+ht);

[min_val, loc_val] = min(x_m);
[max_val, max_loc_val]=max(x_m);
[S_VAL S_LOC] = min(x_e);

max_loc_val=max_loc_val+RLOCATION(k);
loc_val=loc_val+RLOCATION(k);
S_LOC=S_LOC+RLOCATION(k);

x_loc=[x_loc loc_val];
y_val=[y_val min_val];

x_max_loc=[x_max_loc max_loc_val];
y_max_val=[y_max_val max_val];

S_PEAK_LOC=[S_PEAK_LOC S_LOC];
S_PEAK_VAL=[S_PEAK_VAL S_VAL];

end
qq=(S_PEAK_LOC(2)-S_PEAK_LOC(1))/2
    
for k =1:(length(RLOCATION)-1)
x_t=xecg(S_PEAK_LOC(k):RLOCATION(k+1)-qq)

[T_VAL T_LOC] = max(x_t);

T_LOC=T_LOC+S_PEAK_LOC(k);

T_PEAK_LOC=[T_PEAK_LOC T_LOC];
T_PEAK_VAL=[T_PEAK_VAL T_VAL];
end
MM=[];
for k =1:(length(RLOCATION)-1)
    MM(k)=(T_PEAK_LOC(k)+RLOCATION(k+1))/2;
end
for k =1:(length(RLOCATION)-1)
 x_t=xecg(MM(k):RLOCATION(k+1))
 [Q_VAL Q_LOC] = min(x_t);
  Q_LOC=Q_LOC+MM(k);

Q_PEAK_LOC=[Q_PEAK_LOC Q_LOC];
Q_PEAK_VAL=[Q_PEAK_VAL Q_VAL];
end

for k =1:(length(RLOCATION)-1)
   
 x_p=xecg(MM(k):Q_PEAK_LOC(k))
 
[P_VAL P_LOC] = max(x_p);
P_LOC=P_LOC+MM(k);

 P_PEAK_LOC=[P_PEAK_LOC P_LOC];
 P_PEAK_VAL=[P_PEAK_VAL P_VAL];
end

SmidT_VA=[];
SmidT_LO=[];


for k =1:(length(RLOCATION)-1)
   SmidT_L=(S_PEAK_LOC(k)+T_PEAK_LOC(k))/2;
   SmidT_LO(k)= round(SmidT_L)
   SmidT_VA(k) =xecg(SmidT_LO(k))
end



POST_T_VAL=[];
POST_T_LOC=[];

% SmidT_LO=[];
 PLOCATION=[];
 PVALUE=[];
 QLOCATION=[];
 QVALUE=[];
 RRLOCATION=[];
 RRVALUE=[];
 SLOCATION=[];
 SVALUE=[];
 TLOCATION=[];
 TVALUE=[];
 SmidT_VAL=[];
 SmidT_LOC=[];

for i=2:length(RLOCATION)-1
    RRLOCATION(i-1)=RLOCATION(i);
      RRVALUE(i-1)=RVALUE(i);
      SLOCATION(i-1)=S_PEAK_LOC(i);
 SVALUE(i-1)=S_PEAK_VAL(i);
 TLOCATION(i-1)=T_PEAK_LOC(i);
 TVALUE(i-1)=T_PEAK_VAL(i);
   SmidT_LOC(i-1)=SmidT_LO(i);
SmidT_VAL(i-1)=SmidT_VA(i);

end
for i=1:length(RLOCATION)-2
    PLOCATION(i)=P_PEAK_LOC(i);
 PVALUE(i)=P_PEAK_VAL(i);
 QLOCATION(i)=Q_PEAK_LOC(i);
 QVALUE(i)=Q_PEAK_VAL(i);
end   
SmidT_VA=[];
SmidT_LO=[];


for k =1:(length(RLOCATION)-1)
   SmidT_L=(S_PEAK_LOC(k)+T_PEAK_LOC(k))/2;
   SmidT_LO(k)= round(SmidT_L)
   SmidT_VA(k) =xecg(SmidT_LO(k))
end



POST_T_VAL=[];
POST_T_LOC=[];

% SmidT_LO=[];
 PLOCATION=[];
 PVALUE=[];
 QLOCATION=[];
 QVALUE=[];
 RRLOCATION=[];
 RRVALUE=[];
 SLOCATION=[];
 SVALUE=[];
 TLOCATION=[];
 TVALUE=[];
 SmidT_VAL=[];
 SmidT_LOC=[];

for i=2:length(RLOCATION)-1
    RRLOCATION(i-1)=RLOCATION(i);
      RRVALUE(i-1)=RVALUE(i);
      SLOCATION(i-1)=S_PEAK_LOC(i);
 SVALUE(i-1)=S_PEAK_VAL(i);
 TLOCATION(i-1)=T_PEAK_LOC(i);
 TVALUE(i-1)=T_PEAK_VAL(i);
   SmidT_LOC(i-1)=SmidT_LO(i);
SmidT_VAL(i-1)=SmidT_VA(i);

end
for i=1:length(RLOCATION)-2
    PLOCATION(i)=P_PEAK_LOC(i);
 PVALUE(i)=P_PEAK_VAL(i);
 QLOCATION(i)=Q_PEAK_LOC(i);
 QVALUE(i)=Q_PEAK_VAL(i);
end         

POST_T_VA=[];
POST_T_LO=[];
for k =1:(length(RRLOCATION)-1)
  x_t=xecg(TLOCATION(k):PLOCATION(k+1))
  [VAL LOC] = min(x_t);
  LOC=LOC+TLOCATION(k);
POST_T_LO=[POST_T_LO LOC];
POST_T_VA=[POST_T_VA VAL];
end


POST_T_VAL=[];
POST_T_LOC=[];
for k =1:(length(RRLOCATION)-1)
  x_x=(TLOCATION(k)+POST_T_LO(k))/2;
    POST_T_LOC(k)= round(x_x)
   POST_T_VAL(k) =xecg(POST_T_LOC(k))
end




PPG_FOOT_VAL=[];
PPG_FOOT_LOC=[];
PPG_PEAK_LOC=[];
PPG_PEAK_VAL=[];

for i=2:length(RRVALUE)+1
PPG_FOOT_VAL(i-1)=y_val(i);
PPG_FOOT_LOC(i-1)=x_loc(i);
PPG_PEAK_LOC(i-1)=x_max_loc(i);
PPG_PEAK_VAL(i-1)=y_max_val(i);
    
end
dicrotic_val=[];
dicrotic_loc=[];
di_loc=[];
di_val=[];
slope=[];
slop=[];
location=[];
values=[];
y_axis1=[];
y_axis2=[];
for i=1:length(PPG_FOOT_LOC)-1 
    mm=(PPG_FOOT_LOC(i+1)-PPG_PEAK_LOC(i))/2;
  x_f=x2p(PPG_PEAK_LOC(i)+(mm/4):PPG_FOOT_LOC(i+1)-mm),
 %   index = find( abs(diff(x_f)) == 0 )
  y_axis2=x_f(3:length(x_f));
  y_axis1=x_f(1:length(x_f)-2);
 
 slope=y_axis2-y_axis1;
  
  [val pos]=min(abs(slope));
  poss=round(pos)
  value=x_f(poss);
  
  poss=poss+PPG_PEAK_LOC(i)+(mm/4);
   
    values=[values value]
    location=[location poss]
end

% [val dicrotic_loc] =min(slope)
% dicrotic_val=x_f(dicrotic_loc)   
% end
%     
    
    
    
    
    
    
  if select==1
    figure(14)
 subplot(2,1,1)
  plot (t,x, t(RRLOCATION) ,RRVALUE , 'r^',t(SLOCATION),SVALUE,'g*',t(TLOCATION),TVALUE,'r*',QLOCATION/fs,QVALUE,'o',PLOCATION/fs,PVALUE,'O');
 legend('ECG','R peak','S peak' );title('ECG peak detection')
 xlabel('second');ylabel('Volts');
 xlim([0,10]);
 subplot(2,1,2)
  plot([0:length(x2p)-1]/fs,x2p,PPG_FOOT_LOC/fs,PPG_FOOT_VAL,'*r',PPG_PEAK_LOC/fs,PPG_PEAK_VAL,'*g');
 xlabel('second');ylabel('Volts');title('PPG peak detection')
 legend(' PPG finger','PPG FOOT');
 xlim([0 ,10]);
else
    subplot(2,1,1)
  plot ([0:length(x2)-1]/fs,x2,t(RRLOCATION) ,RRVALUE , 'r^',t(SLOCATION),SVALUE,'g*',t(TLOCATION),TVALUE,'r*',QLOCATION/fs,QVALUE,'o',PLOCATION/fs,PVALUE,'O');
 legend('ECG','R peak','S peak' );title('ECG peak detection')
 xlabel('second');ylabel('Volts');
 xlim([0,10]);
 subplot(2,1,2)
  plot([0:length(x2p)-1]/fs,x2p,PPG_PEAK_LOC/fs,PPG_PEAK_VAL,'*r',PPG_FOOT_LOC/fs,PPG_FOOT_VAL,'*g');
 xlabel('second');ylabel('Volts');title('PPG peak detection')
 legend(' PPG finger','PPG FOOT');
 xlim([0 ,10]);
  end

  
   
  if select==1
    figure(24)
 subplot(2,1,1)
  plot (t,x, RRLOCATION/fs ,RRVALUE , 'r^',SLOCATION/fs,SVALUE,'g*',TLOCATION/fs,TVALUE,'r*',QLOCATION/fs,QVALUE,'o',PLOCATION/fs,PVALUE,'O',SmidT_LOC/fs,SmidT_VAL,'*y', POST_T_LOC/fs, POST_T_VAL,'*k');
 legend('ECG','R peak','S peak' );title('ECG peak detection')
 xlabel('second');ylabel('Volts');
 xlim([0,10]);
 subplot(2,1,2)
  plot([0:length(x2p)-1]/fs,x2p,PPG_FOOT_LOC/fs,PPG_FOOT_VAL,'*r',PPG_PEAK_LOC/fs,PPG_PEAK_VAL,'*g',location/fs,values,'*r');
 xlabel('second');ylabel('Volts');title('PPG peak detection')
 legend(' PPG finger','PPG FOOT');
 xlim([0 ,10]);

else
    subplot(2,1,1)
  plot ([0:length(x2)-1]/fs,x2,t(RRLOCATION) ,RRVALUE , 'r^',t(SLOCATION),SVALUE,'g*',t(TLOCATION),TVALUE,'r*',QLOCATION/fs,QVALUE,'o',PLOCATION/fs,PVALUE,'O',SmidT_LOC/fs,SmidT_VAL,'*y', POST_T_LOC/fs, POST_T_VAL,'*k');
 legend('ECG','R peak','S peak' );title('ECG peak detection')
 xlabel('second');ylabel('Volts');
 xlim([0,10]);
 subplot(2,1,2)
  plot([0:length(x2p)-1]/fs,x2p,PPG_PEAK_LOC/fs,PPG_PEAK_VAL,'*r',PPG_FOOT_LOC/fs,PPG_FOOT_VAL,'*g',location/fs,values,'*r');
 xlabel('second');ylabel('Volts');title('PPG peak detection')
 legend(' PPG finger','PPG FOOT');
 xlim([0 ,10]);
  end
PTT_FINGER_PEAK=[];
PTT_FINGER_FOOT=[];
PAT_FINGER_PEAK=[];
PAT_FINGER_FOOT=[];
HR=[]; 
PEP=[];
LVET=[];
EMS=[];
P_AMP=[];
Q_AMP=[];
R_AMP=[];
S_AMP=[];
T_AMP=[];
Q_S_DIST=[];
P_T_DISTANCE=[];
  
  patfingerpeak=[];
   
  patfingerfoot=[];
 
  pttfingerpeak=[];
   
  pttfingerfoot=[];
  sys_area=[];
  dia_area=[];
  total_area=[];
  ratio_area=[];
  pl_ht=[];
  cr_ti=[];
  de_ti=[];
  AI=[];
  RI=[];
  c1=0;
  c2=0;
  c3=0;
  c4=0;
TSN=PPG_PEAK_LOC/fs;
ASN= PPG_PEAK_VAL;
TVN=PPG_FOOT_LOC/fs;
AVN= PPG_FOOT_VAL;
 TDN=  location /fs;
 ADN= values;
 
 
 for i=1:length(ADN)
      x_c=x2p(PPG_PEAK_LOC(i):location(i))
      x_c_up=x_c+1;
  sys_area(i)=sum(x_c_up)
 end
  for i=1:length(ADN)
       x_d=x2p(location(i):PPG_FOOT_LOC(i+1));
       x_d_up=x_d+1;
 dia_area(i)=sum(x_d_up);
 
  end
 for i=1:length(dia_area)
      total_area(i)=sys_area(i)+dia_area(i);
      ratio_area(i)=dia_area(i)/  sys_area(i);
      pl_ht(i)=ASN(i)-AVN(i);
      cr_ti(i)=TSN(i)-TVN(i);
      de_ti(i)=TDN(i)-TSN(i);
      AI(i)=(ADN(i)-AVN(i))/(ASN(i)-AVN(i));
      RI(i)=1-AI(i);
 end
  
for i=1:length(PPG_PEAK_LOC)
  
    patfingerfoot(i)=(PPG_FOOT_LOC(i)/fs)-(RRLOCATION(i)/fs);
   
    patfingerpeak(i)=(PPG_PEAK_LOC(i)/fs)-(RRLOCATION(i)/fs);
    
    if patfingerfoot(i)>0.1 && patfingerfoot(i)<0.7
          PAT_FINGER_FOOT(i)=patfingerfoot(i);
          c1=c1+1;
    else 
        PAT_FINGER_FOOT(i)=0;
    end
        
    if patfingerpeak(i)>0.1 && patfingerpeak(i)<0.7
          PAT_FINGER_PEAK(i)=patfingerpeak(i);
          c2=c2+1;
    else 
        PAT_FINGER_PEAK(i)=0;
    end
    
       
end

PAT_FINGER_FOOT_AVG=sum(PAT_FINGER_FOOT)/c1;
PAT_FINGER_PEAK_AVG=sum(PAT_FINGER_PEAK)/c2;


FINGER_PPG_PEAK_DIST=[];
FINGER_PPG_FOOT_DIST=[];
WRIST_PPG_PEAK_DIST=[];
WRIST_PPG_FOOT_DIST=[];

for i=1:c1-1
    if PAT_FINGER_FOOT(i)>0 & PAT_FINGER_FOOT(i+1)>0
     FINGER_PPG_FOOT_DIST(i)= (PPG_FOOT_LOC(i+1)/fs)-(PPG_FOOT_LOC(i)/fs);
    else
         FINGER_PPG_FOOT_DIST(i)= 0;
    end
end
    
    for i=1:c2-1
    if PAT_FINGER_PEAK(i)>0 & PAT_FINGER_PEAK(i+1)>0
     FINGER_PPG_PEAK_DIST(i)= (PPG_PEAK_LOC(i+1)/fs)-(PPG_PEAK_LOC(i)/fs);
    else
         FINGER_PPG_PEAK_DIST(i)= 0;
    end
    end
      
FINGER_PPG_FOOT_DIST_AVG=sum(FINGER_PPG_FOOT_DIST)/c1;  
FINGER_PPG_PEAK_DIST_AVG=sum(FINGER_PPG_PEAK_DIST)/c2;
PEPP=[];
LVETT=[];
EMSS=[];
for i=1:length(POST_T_LOC)
PEPP(i)=(SmidT_LOC(i)-QLOCATION(i))/fs;
LVETT(i)=(POST_T_LOC(i)-SmidT_LOC(i))/fs;
EMSS(i)=(PEPP(i)+LVETT(i));
end
PEP=sum(PEPP)/length(PEPP);
LVET=sum(LVETT)/length(LVETT);
EMS=sum(EMSS)/length(EMSS);



PTT_FINGER_PEAK_AVG=(sum(PAT_FINGER_PEAK)-(c2*PEP))/c2;
PTT_FINGER_FOOT_AVG=(sum(PAT_FINGER_FOOT)-(c1*PEP))/c1;

HR=7.5*length(RRLOCATION);

a1=abs(PPG_FOOT_VAL(1))+abs(PPG_PEAK_VAL(1));
a2=abs(PPG_FOOT_VAL(2))+abs(PPG_PEAK_VAL(2))
PPG_AMP_RATIO_FINGER=a2/a1;


P_AMP=PVALUE(1);
Q_AMP=QVALUE(1);
R_AMP=RRVALUE(1);
S_AMP=SVALUE(1);
T_AMP=TVALUE(1);

Q_S_DIST=(SLOCATION(1)/fs)-(QLOCATION(1)/fs);
P_T_DISTANCE=(TLOCATION(1)/fs)-(PLOCATION(1)/fs);


prompt={'READING NO.','Enter your SBP:','Enter your DBP:','Enter your REAL HR:'};
% Create all your text fields with the questions specified by the variable prompt.
title='DETAILS'; 
% The main title of your input dialog interface.
answer=inputdlg(prompt,title);
READING_NO = str2num(answer{1});
SB = str2num(answer{2}); 
DB = str2num(answer{3});
HEART = str2num(answer{4});
% Convert these values to a number using str2num.


headings=[];
headings(1)=READING_NO
headings(2)=SB
headings(3)=DB
headings(4)=HEART
headings(5)=PTT_FINGER_PEAK_AVG
headings(6)=PTT_FINGER_FOOT_AVG
headings(7)=PAT_FINGER_PEAK_AVG
headings(8)=PAT_FINGER_FOOT_AVG
headings(9)=HR 
headings(10)=PEP
headings(11)=P_AMP
 headings(12)=Q_AMP
 headings(13)=R_AMP
 headings(14)=S_AMP
 headings(15)=T_AMP
 headings(16)=Q_S_DIST
 headings(17)=P_T_DISTANCE
 headings(18)=PPG_AMP_RATIO_FINGER
 headings(19)=FINGER_PPG_PEAK_DIST(1);
 headings(20)=FINGER_PPG_FOOT_DIST(1)
 headings(21)=FINGER_PPG_PEAK_DIST_AVG
 headings(22)=FINGER_PPG_FOOT_DIST_AVG;
 headings(23)=sys_area(1)
 headings(24)=dia_area(1)
 headings(25)=total_area(1)
 headings(26)=ratio_area(1)
 headings(27)=pl_ht(1)
 headings(28)=cr_ti(1)
 headings(29)=de_ti(1)
 headings(30)=AI(1)
 headings(31)=RI(1)
 headings(32)=EMS
 headings(33)=LVET
 headings(34)=AVN(1)
headings(35)=ASN(1)
 headings(36)=ADN(1)
 
dlmwrite('C.txt',headings,'-append');
