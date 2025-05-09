function harmonizeM()
%translate transmission events from WGS data to colonization probabilities
%output Cpost_M: 7 columns - patient ID, time for colonization probability,
%estimated colonization probability, transmission source patient ID,
%transmission destination patient ID, network distance from the patient to
%source patient.
%load observation data
load('syntheticdata.mat','Y','M','beta_truth')
%Y: clinical culture results. Three columns - sample collection time,
%patient ID, test result (1 positive; 0 negative)
%M: transmission matrix. Mij=1 if WGS data suggest transmission from i to j
%beta_truth: transmission rate used to generate the synthetic outbreak

%load the time range of hospitalization for each patient
load admdis
%admdis: the earliest admission date and last discharge date for each
%patient. Three columns - patient ID, earliest admission data, last
%discharge date. Note each patient may have multiple hospitalizations

%load information on all contacts for each patient during the study period
load neighborinfo
%nlall and partall jointly specify the information of contacts for each
%patient. nlall has three columns - neighbor ID, contact start date,
%contact end date. partall records the rows in nlall corresponding to each
%patient. For instance, the contacts of patient i is
%nlall(partall(i):partall(i+1)-1,1)

%load hospitalizaton duration for each patient
load hospduration
%hospduration: the period of hospitalization for each patient. Three
%columns - patient ID, admission date, discharge date. Note each row is for
%one single hospitalization and one patient may have multiple rows for
%several hospitalizations.

beta=beta_truth;%set the transmission rate used in the inference

N=length(admdis);%total number of patients
Ny=size(Y,1);%number of tested positive

tested=zeros(N,2);
tested(Y(:,2),1)=1;%use this vector to decide which node is tested
tested(Y(:,2),2)=(1:Ny)';%match the rows in M

%set the maximum path lenght L
L=3;
%maximum time
T=50*7+1;

%%%%%%%%%%%%%%%%%%%posterior colonization probability from M
Cpost_M=[];%patient ID, time, posterior, pat1, pat2, distance to pat 1
Cprior=0.5;
Sprior=1-Cprior;

%loop through all pairs
for i=1:Ny
    
    pat1=Y(i,2);
    
    %list of the decendants of pat1
    descendants=unique(Y(M(i,:)>0,2));
    isdes=zeros(N,1);
    isdes(descendants)=1;
    
    if ~isempty(descendants)
        
        [i, Ny]
        tic
        %prepare BFS
        queue=zeros(1e5,1);%create queue
        queue_start=zeros(1e5,1);%the start time of node in paths
        queue_end=zeros(1e5,1);%the end time of node in paths
        front=0;
        tail=0;
        dist=zeros(1e5,1);%distance to pat1 from nodes on paths
        checkedtime=zeros(N,1);%time of first encounter in BFS
        checkeddist=-1*ones(N,1);%distance of first encounter in BFS
        maxdist=0;%current max distance

        %record all paths, auxillary variables
        paths=zeros(1e6,L+1);%record all paths with length <=L
        paths_time_start=zeros(1e6,L+1);%record time for nodes on paths (start time)
        paths_time_end=zeros(1e6,L+1);%record time for nodes on paths (end time)
        paths_contactlength=zeros(1e6,L+1);%record contact length between a node its parent
        pathrows=cell(1e5,1);%match queue, record the rows where node appears
        maxrow=0;%the current max row number in paths
        %check if to keep the paths (containing sequenced patients)
        keeppaths=zeros(1e6,1);
                
        %record unique paths
        uniquepaths=[];%paths with unique patients (paitents may have multiple contact with a person during stay)
        uniquepaths_starttime=[];%the time for first contact in a paths
        uniquepaths_endtime=[];%the time for last contact in a paths (for inference)
        uniquepaths_contactlength=[];%total contact length for a pair of contact

        %%%%%%%%%%%%%%%start BFS
        %push the first element
        tstart=admdis(pat1,2);%first admission of pat1
        checkedtime(pat1)=tstart;
        checkeddist(pat1)=0;
        tail=tail+1;
        queue(tail)=pat1;
        queue_start(tail)=tstart;
        queue_end(tail)=min(admdis(pat1,3),T);
        dist(tail)=0;

        %update paths
        maxrow=maxrow+1;
        paths(maxrow,dist(tail)+1)=pat1;
        paths_time_start(maxrow,dist(tail)+1)=tstart;%first admission of pat1
        paths_time_end(maxrow,dist(tail)+1)=tstart;%include any contacts after admission for seed
        paths_contactlength(maxrow,dist(tail)+1)=0;%no contact with parent as there is no parent
        pathrows{tail}=[pathrows{tail};maxrow];

        while (front~=tail)&&(maxdist<=L)
            %pop the first element
            front=front+1;
            cntpat=queue(front);%current patient
            cntdist=dist(front);%distance of current patient
            
            %find the neighbors of cntpat from checkedtime(cntpat) to T
            %checkedtime(cntpat): first day of cntpat appearing in BFS
            nl_cnt=nlall(partall(cntpat):partall(cntpat+1)-1,:);
            %filter neighbors based on contact time
            if cntpat==pat1 %if the root node, consider all-time neighbors
                %filter neighbors between checkedtime(cntpat) and last discharge
                tstart=checkedtime(cntpat);%the first time of cntpat appearing in paths
                nl_cnt(nl_cnt(:,2)>min(admdis(cntpat,3),T),:)=[];%remove contacts starting after last discharge
                nl_cnt(nl_cnt(:,3)<tstart,:)=[];%remove contacts finished before tstart
            else
                %filter neighbors between the last contact of parent and its discharge in that visit
                tstart=queue_end(front);%last contact of parent
                %check next discharge
                temp=hospduration(hospduration(:,1)==cntpat,:);
                temp=temp((temp(:,2)<=tstart)&(temp(:,3)>=tstart),:);
                if ~isempty(temp)
                    tend=temp(1,3);
                    nl_cnt(nl_cnt(:,2)>tend,:)=[];%remove contacts starting after discharge of current visit
                    nl_cnt(nl_cnt(:,3)<tstart,:)=[];%remove contacts finished before tstart
                else
                    nl_cnt=[];
                end
            end
            
            if ~isempty(nl_cnt)
                neighbors=nl_cnt(:,1);
            else
                neighbors=[];
            end
            
            newencounter=0;
            for j=1:length(neighbors)
                neighbor=neighbors(j);
                tcontact_start=max(nl_cnt(j,2),tstart);%the starting time of contact
                tcontact_end=min(nl_cnt(j,3),T);%the ending time of contact
                if (checkeddist(neighbor)==-1)||(cntdist+1==checkeddist(neighbor))
                    newencounter=newencounter+1;
                end
                
                %first encounter
                if (checkeddist(neighbor)==-1)&&(cntdist+1<=L)
                    %%%%%%push to queue
                    checkedtime(neighbor)=tcontact_start;
                    tail=tail+1;
                    queue(tail)=neighbor;
                    queue_start(tail)=tcontact_start;
                    queue_end(tail)=tcontact_end;
                    dist(tail)=cntdist+1;%include the distance by 1
                    checkeddist(neighbor)=cntdist+1;
                    if dist(tail)>maxdist
                        maxdist=dist(tail);
                    end
                    %%%%%%
                    %update paths, paths_time_start, paths_time_end, paths_contactlength, pathrows
                    pathrows_cnt=pathrows{front};%find the rows of path to update
                    if newencounter==1%the first encountered new neighbor
                        %add neighbor to existing rows
                        for k=1:size(pathrows_cnt,1)
                            %add link only if the new contact end after
                            %the end of last contact
                            if tcontact_end>paths_time_end(pathrows_cnt(k),dist(tail))
                                paths(pathrows_cnt(k),dist(tail)+1)=neighbor;
                                paths_time_start(pathrows_cnt(k),dist(tail)+1)=tcontact_start;
                                paths_time_end(pathrows_cnt(k),dist(tail)+1)=tcontact_end;
                                paths_contactlength(pathrows_cnt(k),dist(tail)+1)=tcontact_end-tcontact_start+1;
                                pathrows{tail}=[pathrows{tail};pathrows_cnt(k)];
                            end
                        end
                        
                    else%addtional neighbors
                        %add new rows in paths
                        for k=1:size(pathrows_cnt,1)
                            if tcontact_end>paths_time_end(pathrows_cnt(k),dist(tail))
                                maxrow=maxrow+1;
                                paths(maxrow,:)=paths(pathrows_cnt(k),:);
                                paths(maxrow,dist(tail)+1)=neighbor;
                                paths_time_start(maxrow,:)=paths_time_start(pathrows_cnt(k),:);%copy previous time
                                paths_time_start(maxrow,dist(tail)+1)=tcontact_start;
                                paths_time_end(maxrow,:)=paths_time_end(pathrows_cnt(k),:);%copy previous time
                                paths_time_end(maxrow,dist(tail)+1)=tcontact_end;
                                paths_contactlength(maxrow,:)=paths_contactlength(pathrows_cnt(k),:);%copy previous time
                                paths_contactlength(maxrow,dist(tail)+1)=tcontact_end-tcontact_start+1;
                                %update pathrows
                                pathrows{tail}=[pathrows{tail};maxrow];
                            end
                            
                        end
                    end
                    %mark the paths to keep (containing sequenced patients)
                    if isdes(neighbor)==1
%                         neighbor
                        keeppaths(pathrows{tail})=1;
                    end
                    
                %first encounter on the same day
                elseif (cntdist+1==checkeddist(neighbor))&&(cntdist+1<=L)
                    %%%%%%%%%%need to update, don't push into queue
                    
                    %%%%%%find neighbor in queue
                    idx=find(queue==neighbor);
                    %%%%%%%%%%%%%%%%
                    %update paths, paths_time, paths_contactlength, pathrows
                    pathrows_cnt=pathrows{front};%find the rows of path to update
                    if newencounter==1%the first encountered new neighbor
                        %add neighbor to existing rows
                        for k=1:size(pathrows_cnt,1)
                            if tcontact_end>paths_time_end(pathrows_cnt(k),dist(idx))
                                paths(pathrows_cnt(k),dist(idx)+1)=neighbor;
                                paths_time_start(pathrows_cnt(k),dist(idx)+1)=tcontact_start;
                                paths_time_end(pathrows_cnt(k),dist(idx)+1)=tcontact_end;
                                paths_contactlength(pathrows_cnt(k),dist(idx)+1)=tcontact_end-tcontact_start+1;
                                pathrows{idx}=[pathrows{idx};pathrows_cnt(k)];
                            end
                        end
                        
                    else%addtional neighbors
                        %add new rows in paths
                        for k=1:size(pathrows_cnt,1)
                            if tcontact_end>paths_time_end(pathrows_cnt(k),dist(idx))
                                maxrow=maxrow+1;
                                paths(maxrow,:)=paths(pathrows_cnt(k),:);
                                paths(maxrow,dist(idx)+1)=neighbor;
                                paths_time_start(maxrow,:)=paths_time_start(pathrows_cnt(k),:);%copy previous time
                                paths_time_start(maxrow,dist(idx)+1)=tcontact_start;
                                paths_time_end(maxrow,:)=paths_time_end(pathrows_cnt(k),:);%copy previous time
                                paths_time_end(maxrow,dist(idx)+1)=tcontact_end;
                                paths_contactlength(maxrow,:)=paths_contactlength(pathrows_cnt(k),:);%copy previous time
                                paths_contactlength(maxrow,dist(idx)+1)=tcontact_end-tcontact_start+1;
                                %update pathrows
                                pathrows{idx}=[pathrows{idx};maxrow];
                            end
                            
                        end
                    end
                    %mark the paths to keep (containing sequenced patients)
                    if isdes(neighbor)==1
%                         neighbor
                        keeppaths(pathrows{idx})=1;
                    end
                    

                end
                
            end
        end
        
        paths=paths(1:maxrow,:);
        paths_time_start=paths_time_start(1:maxrow,:);
        paths_time_end=paths_time_end(1:maxrow,:);
        
        %only keep those with sequenced patients
        paths=paths(keeppaths>0,:);
        paths_time_start=paths_time_start(keeppaths>0,:);
        paths_time_end=paths_time_end(keeppaths>0,:);
        paths_contactlength=paths_contactlength(keeppaths>0,:);
                
        %consolidate paths with the same patients
        checked=zeros(size(paths,1),1);
        for m=1:size(paths,1)
            if checked(m)==0%no overlap found so far
                cntlink=paths(m,:);
                cntlasttime=paths_time_end(m,:);
                cntstarttime=paths_time_start(m,:);
                cntcontactlength=paths_contactlength(m,:);
                for n=m+1:size(paths,1)
                    diff=sum(abs(cntlink-paths(n,:)));
                    if diff==0%match
                        checked(n)=1;%mark found
                        cntlasttime=max(cntlasttime,paths_time_end(n,:));
                        cntstarttime=min(cntstarttime,paths_time_start(n,:));
                        %check if is the same contact. If new, add contact time
                        for k=1:length(cntlink)
                            if (paths_time_start(m,k)~=paths_time_start(n,k))&&(paths_time_end(m,k)~=paths_time_end(n,k))
                                cntcontactlength(k)=cntcontactlength(k)+paths_contactlength(n,k);
                            end
                        end
                        
                    end
                end
                
                %add to unique paths
                uniquepaths=[uniquepaths;cntlink];
                uniquepaths_starttime=[uniquepaths_starttime;cntstarttime];
                uniquepaths_endtime=[uniquepaths_endtime;cntlasttime];
                uniquepaths_contactlength=[uniquepaths_contactlength;cntcontactlength];
                
            end
        end
        
        %%%%%%%%%%%find paths between pairs of transmission
        for k=1:length(descendants)%loop through all descendants
            pat2=descendants(k);
            %find pat2 test time
            Ttest_pat2=Y(Y(:,2)==pat2,1);
            Ttest_pat2=Ttest_pat2(1);
            ispath=zeros(size(uniquepaths,1),1);%if is a path from pat1 to pat2
            for i1=1:size(uniquepaths,1)
                if sum(uniquepaths(i1,:)==pat2)>0
                    %only keep contacts that start before the test of pat2
                    %(transmission occurred before testing pat2
                    if uniquepaths_starttime(i1,uniquepaths(i1,:)==pat2)<Ttest_pat2
                        ispath(i1)=1;
                    end
                end
            end
            %select paths from pat1 to pat2, and relevant info
            cntpaths=uniquepaths(ispath>0,:);
            cntpaths_starttime=uniquepaths_starttime(ispath>0,:);
            cntpaths_endtime=uniquepaths_endtime(ispath>0,:);
            cntpaths_contactlength=uniquepaths_contactlength(ispath>0,:);
            
            %compute path probability
            pathprob=ones(size(cntpaths,1),1);
            %find patients on paths and time
            pats_path=[];%patient ID, last contact time; record patients on paths
            found=zeros(N,1);%whether each patient has been found so far
            for i1=1:size(cntpaths,1)
                idx=find(cntpaths(i1,:)==pat2);
                for i2=2:idx
                    pathprob(i1)=pathprob(i1)*cntpaths_contactlength(i1,i2)*beta;
                end
                for i2=2:idx-1
                    pat3=cntpaths(i1,i2);
                    if found(pat3)==0%new
                        temp=[pat3,cntpaths_endtime(i1,i2),pat1,pat2,i2-1];%pat3, time, pat1, pat2, distance to pat1
                        pats_path=[pats_path;temp];
                    else
                        pats_path(pats_path(:,1)==pat3,2)=max(pats_path(pats_path(:,1)==pat3,2),cntpaths_endtime(i1,i2));
                        pats_path(pats_path(:,1)==pat3,5)=min(pats_path(pats_path(:,1)==pat3,5),i2-1);%shortest distance
                    end
                end
            end
            
            %compute posterior
            for i1=1:size(pats_path,1)
                pat3=pats_path(i1,1);
                %find paths without pat3: H
                Hidx=find(sum(cntpaths==pat3,2)==0);
                %find pthat with pat3: G
                Gidx=find(sum(cntpaths==pat3,2)>0);
                %likelihood: P(pat1->pat2|pat3=S)
                L_S=prod(1-pathprob(Hidx));
                %likelihood: P(pat1->pat2|pat3=C)
                L_C=prod(1-pathprob(Hidx));
                paths_G=cntpaths(Gidx,:);
                pathprob_G=pathprob(Gidx);
                paths_contactlength_G=cntpaths_contactlength(Gidx,:);
                for i2=1:size(paths_G,1)
                    %find the contact length
                    pat3_contactlength=paths_contactlength_G(i2,paths_G(i2,:)==pat3);
                    
                    L_C=L_C*(1-pathprob_G(i2)/(pat3_contactlength*beta));
                end
                
                L_C=1-L_C;
                L_S=1-L_S;
                
                Post_S=Sprior*L_S;
                Post_C=Cprior*L_C;
                
                %normalize to get posterior C
                Post_C=Post_C/(Post_S+Post_C);
                
                temp=[pat3,pats_path(i1,2),Post_C,pats_path(i1,3:5)];%pat3, time, posterior, pat1, pat2, distance
                
                Cpost_M=[Cpost_M;temp];
                
            end
            
        end
        
        toc
    end
    
end

Cpost_M(isnan(Cpost_M(:,3)),:)=[];%remove nan values

%remove duplicated patients, only keep the last estimate (to support more
%inference)
Cpost_M_unique=[];
patchecked=zeros(N,1);
for i=1:size(Cpost_M,1)
    pat=Cpost_M(i,1);
    if patchecked(pat)==0 %not yet encountered
        patchecked(pat)=1;
        temp=Cpost_M(Cpost_M(:,1)==pat,:);
        %find the last estimate
        temp=temp(temp(:,2)==max(temp(:,2)),:);
        %take the average if there are more than one rows (average the
        %evidence)
        temp(1,3)=mean(temp(:,3));
        temp=temp(1,:);%pick the first row
        Cpost_M_unique=[Cpost_M_unique;temp];
    end
end

Cpost_M=Cpost_M_unique;

%pat3, time, posterior, pat1, pat2, distance to pat1
save Cpost_M.mat Cpost_M 