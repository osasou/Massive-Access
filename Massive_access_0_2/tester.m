addpath("utils")
r = 0
re = 0

EbN0 = 100
trials = 100

prop=1
mvalues=[3]
pvalues=[0]
%disp(mvalues)

%%add
patches=2^r;

%tester‚Å‚Í•Ï” K : number of messages‚ğ1‚©‚ç‘‚â‚µ‚Ä‚¢‚éB

%%add

for a = 1:size(mvalues, 2) % a : 1 to m
    %{
    disp('a=')
    disp(a)
    %}
    m=mvalues(a);
    
    for b=1:size(pvalues,2)% b:1to5
        p=pvalues(b);
        i=1;
        %{
        disp('m')
        disp(m)
        disp('p')
        disp(p)
        %}
        %%add
        if (re==0)
            B = m*(m+3)/2 + p;
        else
            B = m*(m+1)/2 + p;
        end
        if (r<1)
            B = patches*B;
        else
            B = patches*B - sum(l(2:end));
        end
        %%add
        output=[];
        K=[];
        dist=[];
        prop=1;

        
        while prop >= 0.0
            K=[K, i];
            disp(i)
            
            [prop, ave_dist] = run(re, m, p, EbN0, i, trials);
            disp("prop")
            disp(prop)
            output = [output, prop];
            dist = [dist, ave_dist];
            
            i=i+1;
            %%add
            
            if(i>5)
                break
            end
            %add
        end
        %output=[output, 0.05];
        %K=[K, i];
        
        filename=strcat("tests/B", num2str(B),'_re', num2str(re),'_m', num2str(m),'_EbN0', num2str(EbN0),'dB_trials', num2str(trials))
        save(filename, "K", "output", "dist");
    end
end



x=strcat("tests/B", num2str(B),'r', num2str(r),'r', num2str(r),'m', num2str(m),'p', num2str(p), 'trials', num2str(trials))
x(1)
plot(K, output)
xlabel('K')
ylabel('output')