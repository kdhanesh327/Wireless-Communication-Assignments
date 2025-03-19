% Dhanesh Kumar

function MQAM(M)

                N=50000;% number of input sumbols
                snrdB=0:0.5:25;
                w=log2(M);
                ENdb=snrdB+log2(w);%ratio of symbol energy and noise in db
                sample=ceil(M.*(rand(N,1)));% generating uniformly distributed 0s and 1s but represented the decimal form of M bits 
                %-----condition for M whether it is perfect square or not------%
                if(ceil(sqrt(M))==floor(sqrt(M)))
                    n=sqrt(M);
                    temp=[];
                    for i=-(n-1):2:n-1  % array formed by both x and y coordinate of constellation will be of same size 
                        for j=-(n-1):2:(n-1)
                            k=i+1i*j;
                            temp(end+1)=k; % M points are stored in array temp
                        end
                    end
                    
                else 
                    temp=[];
                    temp1=[];
                    temp2=[];
                    for i=-(sqrt(2*M)-1):2:(sqrt(2*M)-1)% temp1 is an array of x coordinate of constellation 
                        temp1(end+1)=i;
                    end
                    for j=-(sqrt(0.5*M)-1):2:(sqrt(0.5*M)-1)% temp2 array of y coordinate of constellation
                        temp2(end+1)=j;
                    end
                    for i=1:length(temp1)
                        for j=1:length(temp2)% array of x cordinate and y cordinate have different lengths
                            temp(end+1)=temp1(i)+1i*temp2(j);% constelltion point stored in temp array
                            
                        end
                    end
                   
                end
                
                r=real(temp);
                i=imag(temp);
                Energy=r.^2+i.^2;
                N_Energy=sqrt(sum(Energy)/M);%normalized energy of constellation
                Comp=temp/N_Energy;
                %-----normalized constellation points-----%
                Comp_real=real(Comp);
                Comp_imaginary=imag(Comp);
                hold on
                %---plotting of decision region----%
                
                if (ceil(sqrt(M))==floor(sqrt(M))) % plotting when M is perfect square
                    plot(Comp_real,Comp_imaginary,'o') 
                    xlabel('Inphase')
                    ylabel('Quadrature')
                    title('Constellation Diagram with decision regions')
                    hold on
                    for a=-(n-2)/N_Energy:2/N_Energy:(n-2)/N_Energy
                        xline(a)
                    end
                    for b=-(n-2)/N_Energy:2/N_Energy:(n-2)/N_Energy
                        yline(b)
                    end
                   
                else  %plottting when M is not perfect square
                    figure(1)
                    re=real(temp);
                    im=imag(temp);
                    plot(re/N_Energy,im/N_Energy,'o')
                    xlabel('Inphase')
                    ylabel('Quadrature')
                    title('Constellation Diagram with decision regions')
                    hold on
                    for c=-(sqrt(2*M)-2)/N_Energy:2/N_Energy:(sqrt(2*M)-2)/N_Energy
                        xline(c)
                    end
                    for d=-(sqrt(2*M)-2)/N_Energy:2/N_Energy:(sqrt(2*M)-2)/N_Energy
                        yline(d)
                    end
                   
                end
                 hold off
                
                signal=Comp(sample);%--mapping of msg to MQAM scheme
                signal_i=real(signal);
                signal_q=imag(signal);
                %-----introducing noise------%
                index=1;
               
                for x=ENdb 
                    linearENdb=10.^(x/10);
                    SDnoise=0.707*sqrt(1/(2*linearENdb));%awgn standard noise at different SNR
                    noise=SDnoise*((randn(length(signal),1)+1i*randn(length(signal),1))); %complex noise since our signal is also complex
                    rxd=transpose(signal)+noise;% complex received signal after passing throught awgn channel
                    % minimum distance calculation 
                    rxd_in=real(rxd);
                    rxd_quad=imag(rxd);
                    % repetition of received vector by M times in order to calculate the distancesfrom all the M points
                    rxd_in_rep = repmat(rxd_in,1,M);
                    rxd_quad_rep = repmat(rxd_quad,1,M);
                    distance = zeros(length(rxd_in),M); 
                    min_dist_idx=zeros(length(rxd_in),1);%initializing an array which will store the index of constellation point at which minimum distance occur
                    for l=1:1:N
                        distance(l,:)=(rxd_in_rep(l,:)-Comp_real).^2+ (rxd_quad_rep(l,:)-Comp_imaginary).^2; % calculating the minimum distance
                        [dist,min_dist_idx(l)]=min(distance(l,:));
                    end 
                    y=min_dist_idx; 
                    SEP_calc(index)=vpa((sum(y~=sample)/N),5); % Symbol error probability by doing sum of all those positions where there is difference in message and decoded signal(done through MDD),getting the value upto 10^-5 digit through vpa comman
                    index=index+1;
                    
              end
                EbNolinear=10.^(snrdB/10);
                 % Formulae for symbol error probability
                %SEP_th=1-(1-(1-1/sqrt(M)*erfc(sqrt(3*k*EbNolinear/(2*(M-1)
                %SEP_th=2*(1-(1/sqrt(M)))*erfc(sqrt(1.5*EbNolinear/(M-1)))-((1-(2/sqrt(M))+1/M)*erfc(sqrt(1.5*EbNolinear/(M-1)))*erfc(sqrt(1.5*EbNolinear/(M-1))));
                %SEP_th=1-(1-2*(1-1/sqrt(M)*0.5*erfc(sqrt((3/(M-1))*EbNolinear)))).^2;% calculation of theoritical symbol error probability  
                %SEP_th= 1-(1-(sqrt(M)-1)/sqrt(M))*erfc(sqrt(3/2*k*EbNolinear/(M-1))).^2;
                figure(2)
                plot(snrdB,(SEP_calc),'-')
                hold on
                %---to get SEP_therotical upto 10^-5 value
                z=vpa(1-2*(1-(1/sqrt(M)))*qfunc(sqrt(3*w*EbNolinear/(M-1))),6);
                SEP_theo=vpa((1-z.^2),6);
                %SEP_th=2*(1-(1/sqrt(M)))*erfc(sqrt((3*w*EbNolinear)/(2*(M-1))));
                plot(snrdB,(SEP_theo),'*')
                errorbar(snrdB,(SEP_theo),0.01,"b","LineWidth",1,"LineStyle","-","CapSize",8);
                legend('SEP Calculated','SEP Theoritical')
                xlabel('SNR(dB)');
                ylabel('Symbol Error Rate (Ps)');
                title('SNR (dB) Vs Symbol Error Rate M-QAM scheme');
                grid on;
                hold off
         end
   
       

           