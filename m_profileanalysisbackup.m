Baseelem='manual2';
SEM='Sigma';%Sigma or Quanta

[litfile, litpath] = uigetfile('.txt','MultiSelect','off');
elemdat=readtable("Chemical-Elements-Properties.xlsx");
if strcmp(SEM,'Quanta')
    data=readtable([litpath,litfile],"Delimiter",{'	'},'NumHeaderLines',13,'ReadVariableNames',true);
    vals=data(strcmp(data.InStats_,{'Yes'   }),:);
    err=data(~strcmp(data.InStats_,{'Yes'   }),:);
    varoffset=3;
    endoffset=1;
    signame='Std. deviation';
elseif strcmp(SEM,'Sigma')
    data=readtable([litpath,litfile],'NumHeaderLines',0,'ReadVariableNames',true);
    data.Properties.VariableNames=[data.Properties.VariableNames{1},'Id',data.Properties.VariableNames(2:end-1)];
    vals=data(strcmp(data.Spectrum,{'Spectrum' }),:);
    err=data(~strcmp(data.Spectrum,{'Spectrum' }),:);
    err.Properties.VariableNames=[err.Properties.VariableNames{1},err.Properties.VariableNames(3:end),'NAV'];
    signame={'Sigma'    };
    err=movevars(err,"NAV","After","Spectrum");
    varoffset=3;
    endoffset=1;
end
numelements=numel(vals.Properties.VariableNames);
%% segregation data - if we know which direction things go
manualelem={'Ti','Ta','Ni','Nb','Mo','Fe','Cr','Co','Al','W'};
manualcorr=[-1,1,-1,-1,-1,1,-1,-1,-1,1];
manual=table(manualelem',manualcorr');
manual.Properties.VariableNames={'elem','corr'};
%% DEBUG - plot corrolation plot and a basic plot
%{
figure()
scatter(vals.Ni,vals.Ta)
%}
h=figure()
corrplot(vals(:,3:numelements-1))
figname='CorrPlot';
saveas(h,[litfile(1:end-4),figname,'_base_',Baseelem],'tiffn')
saveas(h,[litfile(1:end-4),figname,'_base_',Baseelem],'fig')

%% rank sort
rank=table();
weightedvalue=rank;
for i=varoffset:numelements-endoffset
    
    elem=vals.Properties.VariableNames{i};
    elemdens=elemdat.Density(strcmp(elemdat.Sym,elem));
    [~,p] = sort(vals.(elem),'ascend');
    r = 1:length(vals.(elem));
    r(p) = r;
    rank.(elem)=r';
    if ~exist("manual")
        linfit=polyfit(vals.(Baseelem),vals.(elem),1);
        if linfit(1,1)>0
            weighteddet=(vals.(elem)-min(vals.(elem)))./table2array(err(strcmp(err.Spectrum,signame),i));%find the weighted interval taking into acound the std in measured values following Ganesan et al 2005
        else 
            weighteddet=(max(vals.(elem))-vals.(elem))./table2array(err(strcmp(err.Spectrum,signame),i));%find the weighted interval taking into acound the std in measured values following Ganesan et al 2005
        end
    else
        if manual.corr(strcmp(manual.elem,elem))==1
            weighteddet=(vals.(elem)-min(vals.(elem)))./table2array(err(strcmp(err.Spectrum,signame),i));%find the weighted interval taking into acound the std in measured values following Ganesan et al 2005
        else 
            weighteddet=(max(vals.(elem))-vals.(elem))./table2array(err(strcmp(err.Spectrum,signame),i));%find the weighted interval taking into acound the std in measured values following Ganesan et al 2005
        end
    end
    weightedvalue.(elem)=weighteddet;
end

averank=mean(table2array(rank),2);%take only the rank data and mean it
aveWeightedRank=mean(table2array(weightedvalue),2);
[~,p] = sort(averank,'ascend');
r = 1:length(averank);
r(p) = r; 
vals.ranksort=r';
[~,p] = sort(aveWeightedRank,'ascend');
r = 1:length(aveWeightedRank);
r(p) = r;
vals.WIRS=r';

%% assign solid fractions
vals.rankFs=(vals.ranksort-0.5)./size(vals,1);
vals.WIRSFs=(vals.WIRS-0.5)./size(vals,1);

%% find density
%first extract density of every element
elemdens=NaN(numelements,1);
for i=varoffset:numelements-endoffset
    elem=vals.Properties.VariableNames{i};
    elemdens(i)=elemdat.Density(strcmp(elemdat.Sym,elem));
end
vals.dens=sum(table2array(vals(:,varoffset:numelements-endoffset))...
    ./100.*repmat(elemdens(varoffset:numelements-endoffset)',size(vals,1),1),2);
%% plot the density as a function of WIRSfs
h=figure();
scatter(vals.WIRSFs,vals.dens)
ylabel('Density (g/cm^3)')
xlabel('Fs')
title('WIRS Density')
figname='WIRS Density';
saveas(h,[litfile(1:end-4),figname,'_base_',Baseelem],'tiffn')
saveas(h,[litfile(1:end-4),figname,'_base_',Baseelem],'fig')
%% plot compositions as a function of Fs
h=figure();
scatter(vals.WIRSFs,table2array(vals(:,3:numelements-1)))
ylabel('Wt %')
xlabel('Fs')
legend(vals.Properties.VariableNames{3:numelements-1})
title('WIRS')
figname='WIRS Composition TOT';
saveas(h,[litfile(1:end-4),figname,'_base_',Baseelem],'tiffn')
saveas(h,[litfile(1:end-4),figname,'_base_',Baseelem],'fig')
%% another type of plot
h=figure();
for i=varoffset:numelements-endoffset
    subplot(3,4,i-varoffset+1)
    scatter(vals.WIRSFs,table2array(vals(:,i)))
    ylabel('Wt %')
    xlabel('Fs')
    title(vals.Properties.VariableNames{i})
end
sgtitle('WIRS')
figname='WIRS';
saveas(h,[litfile(1:end-4),figname,'_base_',Baseelem],'tiffn')
saveas(h,[litfile(1:end-4),figname,'_base_',Baseelem],'fig')
h=figure();
for i=varoffset:numelements-endoffset
    subplot(3,4,i-varoffset+1)
    scatter(vals.rankFs,table2array(vals(:,i)))
    ylabel('Wt %')
    xlabel('Fs')
    title(vals.Properties.VariableNames{i})
end
subplot(3,4,i-varoffset+2)
scatter(vals.rankFs,table2array(vals(:,i)))
ylabel('Wt %')
xlabel('Fs')
title(vals.Properties.VariableNames{i})
sgtitle('Rank Sort')
figname='Rank_Sort';
saveas(h,[litfile(1:end-4),figname,'_base_',Baseelem],'tiffn')
saveas(h,[litfile(1:end-4),figname,'_base_',Baseelem],'fig')