% Dhanesh Kumar

function MPSK(M)
                

                snrdB=0:0.5:25;%range of ration of bit energy by noise
                k=log2(M);
                ENdb=snrdB+10*log10(k);% conversion of bit energy by noise ration to symbol energy by noise ratio
                N=50000;%number of input symbols
                sample=ceil(M.*(rand(N,1)));% generating uniformly distributed 0s and 1s but represented the decimal form of M bits 
                %---mapping of message to mpsk scheme----%
                Inphase=1/sqrt(2)*(cos (2*pi*(sample-1)/M));
                Quadrature =1/sqrt(2)*(sin(2*pi*(sample-1)/M));
                M_psk=Inphase+1i*Quadrature;% M-psk signal
                %----mpsk constellation points-----%
                signal_in=zeros(1,M);
                signal_quad=zeros(1,M);
                for j=1:1:M
                    signal_in(j)=1/sqrt(2)*cos((j-1)/M*2*pi);
                    signal_quad(j)=1/sqrt(2)*sin((j-1)/M*2*pi);
                end
                figure(1)
                plot(signal_in,signal_quad,'s')
                %plotting constellation points
                hold on
                %-----plotting of deccision region------%
                a=linspace(0,2*pi,8192);
                r=1;
                x=r*cos(a+(pi/M));
                y=r*sin(a+(pi/M));
                axis equal;
                plot([zeros(1,M); x(1:(8192/M):end)], [zeros(1,M); y(1:(8192/M):end)])
                xlabel('Inphase')
                ylabel('Quadrature')
                title('Constellation Diagram')
                hold off
                SEP_calc=zeros(1,length(ENdb));
                %-----introducing noise------%
                index=1;
                for x=ENdb
                    linearENdb=10.^(x/10);
                    SDnoise=0.707*sqrt(1/(2*linearENdb));%awgn standard noise at different SNR
                    noise=SDnoise*((randn(length(M_psk),1)+1i*randn(length(M_psk),1)));%complex noise since our signal is also complex
                    rxd=M_psk+noise;% complex received signal after passing throught awgn channel
                    % minimum distance calculation 
                    rxd_i=real(rxd);
                    rxd_q=imag(rxd);
                    % repetition of received vector by M times in order to calculate the distancesfrom all the M points
                    rxd_i_rep = repmat(rxd_i,1,M);
                    rxd_q_rep = repmat(rxd_q,1,M);
                    distance = zeros(length(rxd_i),M); 
                    min_dist_idx=zeros(length(rxd_i),1);%initializing an array which will store the index of constellation point at which minimum distance occur
                    for l=1:1:N
                        distance(l,:)=(rxd_i_rep(l,:)-signal_in).^2+ (rxd_q_rep(l,:)-signal_quad).^2;
                        [dist,min_dist_idx(l)]=min(distance(l,:));
                    end
                    y=min_dist_idx;
                    
                  SEP_calc(index)=vpa((sum(y~=sample)/N),5);% Symbol error probability by doing sum of all those positions where there is difference in message and decoded signal(done through MDD),getting the value upto 10^-5 digit through vpa comman
                  index=index+1;
                    
                end
                EbNolinear=10.^(snrdB/10);
                SEP_theo=vpa((erfc(sqrt(EbNolinear*k)*sin(pi/M))),5);% calculation of theoritical symbol error probability and to get SEP_therotical upto 10^-5 value
                figure(2)
                semilogy(snrdB,(SEP_calc),'-')
                hold on
                semilogy(snrdB,(SEP_theo),'*')
                legend('SEP Calculated','SEP Theoritical')
                xlabel('SNR(dB)');
                ylabel('Symbol Error Rate (Ps)');
                title('SNR (dB) Vs Symbol Error Rate M-PSK scheme');
                grid on;
                hold off
         end
         
    
            