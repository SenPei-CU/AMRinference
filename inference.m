function inference()
%infer patient colonization probability
%output xpost: estimated colonization probability for each patient on each
%day t, recorded in xpost{t}. Each row corresponds to each patient (patient
%ID matching patlist{t}), each column is the estimated colonization
%probability for each ensemble member.
%compile C++ before running the function
% mex ABMsimulation_chain_p_beta.c
% mex masterequation_p_beta.c

%load contact network
load contactnetwork
%patlist: the list of patient IDs on each day
%adm: the first admission date for each patient
%readmittedpat: patient IDs of readmitted patients on each day
%nl and part jointly specify the contact network structure on each day.
%Specifically, nl{t} is the list of contact IDs on day t. part(:,t) records
%the rows of contacts for each patient in nl{t}. The neighbors of patient i
%on day t are nl{t}(part(i):part(i+1)-1,t)

%load observations and parameters
load('syntheticdata.mat','Y','gamma_truth','beta_truth','alpha_truth','rho')
%load colonization probabilities estimated using WGS
load Cpost_M

%prepare variables
num_ens=100;%ensemble member
num_pat=size(adm,1);
T=50*7;%number of days
T_week=floor(T/7);%number of weeks
x=-1*ones(num_pat,num_ens);%colonized probability: -1: before first hospitalization
xpost=cell(T,1);

%set up parameters
gamma_ref=gamma_truth*ones(1,num_ens);%reference importation rate
beta=beta_truth*ones(1,num_ens);%transmission rate
alpha=alpha_truth*ones(num_pat,num_ens);%decolonization rate

%use regression model to initialize model
gamma_pat=zeros(num_pat,num_ens);
%load medicalinfo
load medicalinfo
%only keep first adm in medicalinfo
medicalinfo=medicalinfo(medicalinfo.FirstAdm==1,:);
%load trained model, using data prior to 2012
load regressionmodel_select
%prepare variables (14, full model)
xvar=[];
colname=medicalinfo.Properties.VariableNames;
for i=1:length(variables)
    idx=find(strcmp(colname,variables{i})>0);
    xvar=[xvar,medicalinfo{:,idx}];
end
y = glmval(bs(:,ceil(rand(1,num_ens)*size(bs,2))),xvar,'logit');

%rescale y=y^a to match the mean importation rate
ytemp=reshape(y,[],1);
%find the best a
a=1:0.1:10;
ymean=zeros(length(a),1);
for j=1:length(a)
    ymean(j)=mean(ytemp.^a(j));
end
diff=abs(ymean-gamma_ref(1));
besta=a(diff==min(diff));

for i=1:num_ens
    y(:,i)=y(:,i).^besta;
end

idx=medicalinfo.PatID;
y=y(idx>0,:);
idx=idx(idx>0);
gamma_pat(idx,:)=y;
gamma=gamma_pat;%importation rate
para=[gamma;beta];
%set initial conditions
inipat=patlist{1}(:,1);
x(inipat,:)=gamma(inipat,:);%colonization probability for patients initially in hospital
%start from day 0
%loop through time, update every week
for t=1:T_week%update for 50 weeks
    [t,T_week]
    tic
    tcnt=(t-1)*7+1;%first day of week t
    nl_cnt=nl{tcnt};
    part_cnt=part(:,tcnt);
    
    %find current hospitalized patients
    patlist_cnt=patlist{tcnt}(:,1);
    %find nodes with positive probability
    validnode=find(x(:,1)>0);
    patlist_cnt=intersect(patlist_cnt,validnode);
    
    %find nodes that are tested after tcnt
    pat_tested=intersect(patlist_cnt,Y(Y(:,1)>tcnt,2));
    info_test=Y(ismember(Y(:,2),pat_tested),[2,1,3]);%patID,test time,result
    
    %find nodes on paths
    pat_onpath=intersect(patlist_cnt,Cpost_M(Cpost_M(:,2)>tcnt,1));
    info_path=Cpost_M(ismember(Cpost_M(:,1),pat_onpath),[1:3]);%patID, onpath time, posterior

    %get likelihood for tested
    if ~isempty(info_test)
        [L1_C,L1_S]=getlikelihood_test_fast(tcnt,info_test,x,nl,part,patlist,para,alpha,readmittedpat,rho);
        %check if L1 contains nan
        if sum(isnan(L1_C))+sum(isnan(L1_S))>0
            disp('NaN in L1!');
            break
        end
    end
    %get likelihood for patients on paths
    if ~isempty(info_path)
        [L2_C,L2_S]=getlikelihood_path_fast(tcnt,info_path,x,nl,part,patlist,para,alpha,readmittedpat);
        %check if L2 contains nan
        if sum(isnan(L2_C))+sum(isnan(L2_S))>0
            disp('NaN in L2!');
            break
        end
    end
    
    %update nodes and their neighbors
    pat_update=union(pat_tested,pat_onpath);%nodes need to be updated
    
    for i=1:length(pat_update)
        pat=pat_update(i);
        
        %compute the posterior
        Cprior=x(pat,:);
        
        if sum(isnan(Cprior))>0
            disp('NaN in Cprior!');
            break
        end
        
        %if var(Cprior)==0, add nosie
        if var(Cprior)==0
            Cprior=Cprior+randn(1,num_ens)*0.01;
            Cprior(Cprior>1)=1;%check range
            Cprior(Cprior<0)=0;
        end
        prior_var=max(var(Cprior),1e-4);%make sure prior_var>0
        %likelihood
        LL_C=1;%likelihood for xi(t)=C
        LL_S=1;%likelihood for xi(t)=S
        %check if pat is in test
        if ismember(pat,pat_tested)
            idx=find(pat_tested==pat);
            LL_C=LL_C*L1_C(idx,:);
            LL_S=LL_S*L1_S(idx,:);
        end
        %check if pat is in onpath
        if ismember(pat,pat_onpath)
            idx=find(pat_onpath==pat);
            LL_C=LL_C*L2_C(idx,:);
            LL_S=LL_S*L2_S(idx,:);
        end
        
        if sum(isnan(LL_C))>0
            disp('NaN in LL_C!');
            break
        end
        
        if sum(isnan(LL_S))>0
            disp('NaN in LL_S!');
            break
        end
        
        Cpost=Cprior.*LL_C;
        Spost=(1-Cprior)*LL_S;
        Cpost=Cpost./(Cpost+Spost);
        %check if there are nans
        if sum(isnan(Cpost))>0
            disp('NaN in Cpost!')
            Cpost(isnan(Cpost))=mean(Cpost(~isnan(Cpost)));
        end
        %check range
        Cpost(Cpost>1)=1;
        Cpost(Cpost<0)=0;
        %update the focal node
        dy=Cpost-Cprior;
        x(pat,:)=Cpost;
        %loop through neighbors
        for j=part_cnt(pat):part_cnt(pat+1)-1
            nei=nl_cnt(j);
            nei_prior=x(nei,:);
            r=cov(Cprior,nei_prior);
            r=r(2,1)/prior_var;
            dx=r*dy;
            if sum(isnan(dx))>0
                dx=0;
            end
            
            temp=nei_prior+dx;
            %check range
            temp(temp>1)=1;
            temp(temp<0)=0;
            %update
            x(nei,:)=temp;
            
            if sum(isnan(temp))>0
                disp('NaN in x(nei,:)!');
                break
            end
            
        end
    end
    
    %record x posterior
    xpost{tcnt}=x(patlist_cnt,:);
    
    %integrate to next time step (one week)
    for t1=tcnt:tcnt+6
        nl_cnt=nl{t1};
        part_cnt=part(:,t1);
        patlist_cnt=patlist{t1}(:,1);
        degree=zeros(num_pat,1);
        degree(patlist_cnt)=patlist{t1}(:,2);
        readmittedpat_cnt=readmittedpat{t1};%patients outside hospital who are readmitted later
        %%%%%%%%%run master equations
        
        x=masterequation_p_beta(x,nl_cnt,part_cnt,patlist_cnt,degree,para,alpha,readmittedpat_cnt);

        %record x posterior
        xpost{t1}=x(patlist_cnt,:);
    end
    toc
end

save xpost.mat xpost


function [L1_C,L1_S]=getlikelihood_test_fast(tcnt,info_test,x,nl,part,patlist,para,alpha,readmittedpat,rho)
%rho is the effective sensitivity
n=size(info_test,1);
num_pat=size(part,1)-1;
L1_C=zeros(n,1);%likelihood of observation if current status is C
L1_S=zeros(n,1);%likelihood of observation if current status is S

for i=1:n
    pat=info_test(i,1);
    Ttest=info_test(i,2);%testing date
    %test results show status at beginning of Ttest (end of Ttest-1)
    %%%%%%%%%%%%%%%%%%set status to C at tcnt
    xtemp=mean(x,2);%record current status
    xtemp(pat,:)=1;%set the status to C
    for t=tcnt:Ttest-1
        nl_cnt=nl{t};
        part_cnt=part(:,t);
        patlist_cnt=patlist{t}(:,1);
        degree=zeros(num_pat,1);
        degree(patlist_cnt)=patlist{t}(:,2);
        readmittedpat_cnt=readmittedpat{t};%patients outside hospital who are readmitted later
        %%%%%%%%%run master equations
        xtemp=masterequation_p_beta(xtemp,nl_cnt,part_cnt,patlist_cnt,degree,para,alpha,readmittedpat_cnt);
    end
    %get colonization probability at time Ttest
    PC_C=xtemp(pat,:);
    %%%%%%%%%%%%%%%%%%set status to S at tcnt
    xtemp=mean(x,2);%record current status
    xtemp(pat,:)=0;%set the status to S
    for t=tcnt:Ttest-1
        nl_cnt=nl{t};
        part_cnt=part(:,t);
        patlist_cnt=patlist{t}(:,1);
        degree=zeros(num_pat,1);
        degree(patlist_cnt)=patlist{t}(:,2);
        readmittedpat_cnt=readmittedpat{t};%patients outside hospital who are readmitted later
        %%%%%%%%%run master equations
        xtemp=masterequation_p_beta(xtemp,nl_cnt,part_cnt,patlist_cnt,degree,para,alpha,readmittedpat_cnt);
    end
    %get colonization probability at time Ttest
    PC_S=xtemp(pat,:);
    %%%%%%%%%%%%%%%%%%compute likelihood
    %L1_C
    L1_C(i,:)=rho*PC_C;
    %L1_S
    L1_S(i,:)=rho*PC_S;
end

function [L2_C,L2_S]=getlikelihood_path_fast(tcnt,info_path,x,nl,part,patlist,para,alpha,readmittedpat)
n=size(info_path,1);
num_pat=size(part,1)-1;
L2_C=zeros(n,1);%likelihood of observation if current status is C
L2_S=zeros(n,1);%likelihood of observation if current status is S

for i=1:n
    pat=info_path(i,1);
    Ttest=info_path(i,2);%hospitalization date
    %test results show status at beginning of Ttest (end of Ttest-1)
    %%%%%%%%%%%%%%%%%%set status to C at tcnt
    xtemp=mean(x,2);%record current status
    xtemp(pat,:)=1;%set the status to C
    for t=tcnt:Ttest-1
        nl_cnt=nl{t};
        part_cnt=part(:,t);
        patlist_cnt=patlist{t}(:,1);
        degree=zeros(num_pat,1);
        degree(patlist_cnt)=patlist{t}(:,2);
        readmittedpat_cnt=readmittedpat{t};%patients outside hospital who are readmitted later
        %%%%%%%%%run master equations   
        xtemp=masterequation_p_beta(xtemp,nl_cnt,part_cnt,patlist_cnt,degree,para,alpha,readmittedpat_cnt);
    end
    %get colonization probability at time Ttest
    PC_C=xtemp(pat,:);
    %%%%%%%%%%%%%%%%%%set status to S at tcnt
    xtemp=mean(x,2);%record current status
    xtemp(pat,:)=0;%set the status to S
    for t=tcnt:Ttest-1
        nl_cnt=nl{t};
        part_cnt=part(:,t);
        patlist_cnt=patlist{t}(:,1);
        degree=zeros(num_pat,1);
        degree(patlist_cnt)=patlist{t}(:,2);
        readmittedpat_cnt=readmittedpat{t};%patients outside hospital who are readmitted later
        %%%%%%%%%run master equations
        xtemp=masterequation_p_beta(xtemp,nl_cnt,part_cnt,patlist_cnt,degree,para,alpha,readmittedpat_cnt);
    end
    %get colonization probability at time Ttest
    PC_S=xtemp(pat,:);
    %%%%%%%%%%%%%%%%%%compute likelihood
    %L2_C
    L2_C(i,:)=info_path(i,3)*PC_C+(1-info_path(i,3))*(1-PC_C);
    %L1_S
    L2_S(i,:)=info_path(i,3)*PC_S+(1-info_path(i,3))*(1-PC_S);
end