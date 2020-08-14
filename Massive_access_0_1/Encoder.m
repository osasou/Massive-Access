classdef Encoder
%Class responsible for encoding chirps - contains all functons pertaining to encoding
    properties
        r               % 0, 1 or 2; 2^r patches
        re              % logical: false(0)=complex chirps, true(1)=real chirps
        m               % size of chirp code (2^m is length of codewords in each slot), recommend m = 6, 7 or 8
        p               % 2^p slots require p<=(m-r)(m-r+3)/2-1 for complex
                        % p<=(m-r)(m-r+1)/2-1 for real
        K               % number of messages
        EbN0            % energy-per-bit (Eb/N0)
        input_bits      % raw input bitstring
        B               % number of bits being encoded
        patches         % number of patches
    end

    methods
        function self = Encoder(re,m,p,K,EbN0,input_bits)
            addpath('utils');
            self.r = 0;
            self.re = re;
            self.m = m;
            self.p = p;
            self.K = K;
            self.EbN0 = EbN0;
            self.input_bits = input_bits;
            self.patches=2^self.r;
            if (re==0)
                 self.B = m*(m+3)/2 + p;
            else
                 self.B = m*(m+1)/2 + p;
            end
        end

        function [self,bits] = generate_random_bits(self)
        % generates some random bits to pass into encoder
        % row 行 : B(number of bits being encoded)
        % column 列 : K(number of messages)
            bits = rand(self.B,self.K) > 0.5; 
            %disp("bits")
            %disp(bits)
            self.input_bits=bits;
        end


        function [Y, h_all] = chirrup_encode(self)

        %chirrup_encode  Generates K random messages and performs CHIRRUP encoding
        %
        % Y            Y{p} is a 2^m x 2^p matrix of measurements for patch p
        % input_bits   B x K matrix of the K B-bit messages
        % parity       parity check codes generated in the tree encoding
        %
        % No. of messages is B = 2^r*[(m-r-p)(m-r-p+3)/2+p-1]-sum(l)  for complex
        %                          B = 2^r*[(m-r-p)(m-r-p+1)/2+p-1)-sum(l) for real
        %
        % AJT (12/9/18)

            global h_all
            
            %generate random messages
            patch_bits = self.input_bits.';
            
            %generate measurements for each patch
            for patch = 1:self.patches
                sigma = sqrt(self.patches*2^self.m/(self.B*self.EbN0));
                [Y{patch}, h_all] = self.sim_from_bits(sigma,patch_bits(:,:,patch));
            end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [Y, h_all] = sim_from_bits(self,sigma,bits)

        % sim_from_bits  Generates binary chirp measurements（大きさ） from bits
        % sigma      SD of noise: sigma = sqrt(patches*2^m/(B*EbN0))
        % bits       k x 2^m matrix of bits to encode
        %
        % Y          length-2^m vector of measurements
        %
        % AJT (12/9/18)
        h_all = [];
        alpha = 1.0;
        circle_r = 1;
        active_user_num = self.K;
        
            Y = zeros(2^self.m,2^self.p);
            h_all = [];
            %Yは2^p個あるslotsの各受信語を示している
            
            %p=0なので、zeros(2^m,1)となる
            for k = 1:self.K %the number of active user で回す
                [x,y] = random_circle(circle_r, 1);
                d_dash = abs(x).^2 + abs(y).^2;
                d = sqrt(d_dash);
                d = 1;
                
                bits1 = [bits(k,:)];
%                 bits1 = [1 1 1 1 1 1 0 0 0];
                [Pee1,bee1] = self.makePb(bits1);
                
                %generate binary chirp vector for each slot
                rm1 = self.gen_chirp(Pee1,bee1);

                h = normrnd(0,0.5)+1i*normrnd(0,0.5);
                 h = 1 + 1i;
                  h = 1;
%                 h = 1;
%                 h = 1i;
                h_all=[h_all, h];
                
                Y(:,1) = Y(:,1)+h*rm1*1/(d^alpha);
            end
%                sigma = 0;
            %add noise (Gaussian for real, Complex Gaussian for complex)
            if (self.re==0)
                %B = repmat(A,n) は、行と列の次元に A のコピーを n 個含む配列を返します。
                Y = Y + repmat(sigma*(randn(2^self.m,1)+1i*randn(2^self.m,1)),[1 2^self.p]);
            else
                Y = Y + repmat(sigma*randn(2^self.m,1),[1 2^self.p]);
            end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [P,b] = makePb(self,bits)
        % generates a P and b from a bit string
        %
        % bits      vector of bits
        %
        % P     symettric real matrix
        % b     real vector

            if (self.re==0)
                nMuse = self.m*(self.m+1)/2;
            else
                nMuse = self.m*(self.m-1)/2;
            end
            basis = makeDGC(self.re,self.m);
            Pbits = bits(1:nMuse); %bits[1:3]の前半分
            
            P = mod( sum(basis(:,:,find(Pbits)),3), 2);
            %sum(A,3):行列Aの３番目の次元に沿って和を計算する。つまり、行列の足し算なだけ
            b = bits(nMuse+1:nMuse+self.m);
            % bits[4:6]の後ろ半分
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    methods(Static)

            function rm = gen_chirp(P,b)
            % generates a read-muller code from an input P and b

                M = length(b);
                rm = zeros(2^M,1);
                a = zeros(M,1);
                for q = 1:2^M
                    sum1 = a'*P*a;
                    sum2 = b*a;
                    rm(q) = i^sum1 * (-1)^sum2;
                    % next a
                    for ix = M:-1:1
                        if a(ix)==1
                            a(ix)=0;
                        else
                            a(ix)=1;
                            break;
                        end
                    end
                end
            end
    end


end