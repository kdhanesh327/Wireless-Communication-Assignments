% Dhanesh Kumar

function main2(filename)
 tic
                table1=readtable(filename);% reading data from csv file which contains the arbitrary 4 points of constellation
                temp=table2array(table1);
                snrdb=0:2:20;
                M=4;
                k=log2(M);
                ENdb=snrdb+10*log10(k); %ratio of symbol energy and noise in db
                N=500000;% number of input sumbols
                %-----plotting of decision region ---%
                
                x_c=temp(:,1);
                y_c=temp(:,2);
                %----making a regular plane---%
                pmax=max(x_c)+50;
                pmin=min(x_c)-50;
                qmax=max(y_c)+50;
                qmin=min(y_c)-50;
                %------taking large sample of points of defined plane-----%
                point1=[];
                for i=pmin:1:pmax
                    for j=qmin:1:qmax
                        k=i+1i*j;
                        point1(end+1)=k;
                    end
                end
                point=transpose(point1);
                x_cpoint=real(point);
                y_cpoint=imag(point);
                x_cpoint_rep=repmat(x_cpoint,1,length(x_c));
                y_cpoint_rep=repmat(y_cpoint,1,length(y_c));
                %calculating the distances of all the points taken on a plane
                %with the given constellation points.
                
                distance=zeros(length(pmax),length(x_c));
                min_dist_idx=zeros(length(point),1);
                %--segregating the points according to the minimum distance 
                b1=[];
                b2=[];
                r1=[];
                r2=[];
                g1=[];
                g2=[];
                k1=[];
                k2=[];
                %---assigning 4 different colours to define decision region
                for l=1:1:length(point)
                    distance(l,:)=(x_cpoint_rep(l,:)-transpose(x_c)).^2+(y_cpoint_rep(l,:)-transpose(y_c)).^2;
                    [D,min_dist_idx(l)]=min(distance(l,:));
                    if (min_dist_idx(l)==1)
                        b1(end+1)=x_cpoint(l);
                        b2(end+1)=y_cpoint(l);
                    end
                    if (min_dist_idx(l)==2)
                        r1(end+1)=x_cpoint(l);
                        r2(end+1)=y_cpoint(l);
                    end
                    if (min_dist_idx(l)==3)
                        g1(end+1)=x_cpoint(l);
                        g2(end+1)=y_cpoint(l);
                    end
                    if (min_dist_idx(l)==4)
                        k1(end+1)=x_cpoint(l);
                        k2(end+1)=y_cpoint(l);
                    end
                end
                figure(1)
                plot(b1,b2,'b.')
                hold on
                plot(r1,r2,'r.')
                hold on
                plot(g1,g2,'g.')
                hold on
                plot(k1,k2,'k.') 
                plot(x_c,y_c,'d')
                str={'point1','point2','point3','point4'};
                text(x_c,y_c,str)
                hold off
                xlabel('Inphase')
                ylabel('Quadrature')
                title('Constellation Diagram with decision region')
                legend('point1','point2','point3','point4')
                c_p=x_c+1i*y_c;% constellation points
                c_p_t=c_p';
                msg=ceil(M.*(rand(N,1)));% generating uniformly distributed 0s and 1s but represented the decimal form of M bits 
                signal=c_p_t(msg);%mapping of msg to contellation points of arbitrary scheme
                signal_i=real(signal);
                signal_q=imag(signal);
                SEP_calc=zeros(1,length(ENdb));
                %-----introducing noise------%
                index=1;
                for x=ENdb
                    linearENdb=10.^(x/10);
                    SDnoise=1/sqrt(2)*(sqrt(1/(2*linearENdb)));%awgn standard noise at different SNR
                    noise=5*(SDnoise*((randn(length(signal),1)+1i*randn(length(signal),1))));%complex noise since our signal is also complex
                    rxd=transpose(signal)+noise;% complex received signal after passing throught awgn channel
                    rxd_i=real(rxd);
                    rxd_q=imag(rxd);
                    % repetition of received vector by M times in order to calculate the distancesfrom all the M points
                    rxd_i_rep = repmat(rxd_i,1,M);
                    rxd_q_rep = repmat(rxd_q,1,M);
                    distance = zeros(length(rxd_i),M); 
                    min_dist_idx=zeros(length(rxd_i),1);%initializing an array which will store the index of constellation point at which minimum distance occur
                    for l=1:1:N
                        distance(l,:)=(rxd_i_rep(l,:)-transpose(x_c)).^2+ (rxd_q_rep(l,:)-transpose(y_c)).^2;
                        [D,min_dist_idx(l)]=min(distance(l,:));
                    end
                    y=min_dist_idx;
                    SEP_calc(index)=sum(y~=msg)/N;% Symbol error probability by doing sum of all those positions where there is difference in message and decoded signal(done through MDD)
                    index=index+1;
                end
                EbNolinear=10.^(snrdb/10);
               
                figure(2)
                plot(snrdb,(SEP_calc),'-')
                xlabel('SNR(dB)');
                ylabel('Symbol Error Rate (Ps)');
                title('SNR (dB) Vs Symbol Error Rate arbitrary  scheme');
                grid on;
                hold off;
 toc
 end           