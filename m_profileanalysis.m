%Core variables to adjust
%---------------------------------------------------------------------------------------
Baseelem='Al';
SEM='Tescan';%Sigma or Quanta or Tescan
%----------------------------------------------------------------------------------------------
%{
%% segregation data - if we know which direction things go

manualelem={'Ti','Ta','Ni','Nb','Mo','Fe','Cr','Co','Al','W'};
manualcorr=[-1,1,-1,-1,-1,1,-1,-1,-1,1];
manual=table(manualelem',manualcorr');
manual.Properties.VariableNames={'elem','corr'};
%}

% read the files to create variables
[EDSfile, EDSpath] = uigetfile('.txt','MultiSelect','off');%select EDS file
elemdat=readtable("Chemical-Elements-Properties.xlsx");
%%
if strcmp(SEM,'Quanta')
    %read the file as a table
    data=readtable([EDSpath,EDSfile],"Delimiter",{'	'},'NumHeaderLines',13,'ReadVariableNames',true);
    vals=data(strcmp(data.InStats_,{'Yes'   }),:);
    err=data(~strcmp(data.InStats_,{'Yes'   }),:);
    varoffset=3;
    endoffset=1;
    signame='Std. deviation';
elseif strcmp(SEM,'Sigma')
    data=readtable([EDSpath,EDSfile],'NumHeaderLines',0,'ReadVariableNames',true);
    data.Properties.VariableNames=[data.Properties.VariableNames{1},'Id',data.Properties.VariableNames(2:end-1)];
    vals=data(strcmp(data.Spectrum,{'Spectrum' }),:);
    err=data(~strcmp(data.Spectrum,{'Spectrum' }),:);
    err.Properties.VariableNames=[err.Properties.VariableNames{1},err.Properties.VariableNames(3:end),'NAV'];
    signame={'SigmaMean'};
    err=movevars(err,"NAV","After","Spectrum");
    varoffset=3;
    endoffset=1;
elseif strcmp(SEM,'Tescan')
    data=readtable([EDSpath,EDSfile],'NumHeaderLines',7,'ReadVariableNames',true);
    data2 = table2array(data(:,2:end));
    data2 = array2table(data2.');
    data2.Properties.VariableNames = data.SpectrumLabel;
    vals=data2(:,1:end-1);
    [a,b]=find(isnan(vals{:,:}));
    vals{a,b}=0;
    varoffset=1;
    endoffset=0;
    err=readtable([EDSpath,EDSfile],'NumHeaderLines',1,'ReadVariableNames',true);
    err=err(1:4,:);
    err2 = table2array(err(:,2:6));
    err2 = array2table(err2.');
    err2.Properties.VariableNames = err.Statistics;
%     err2.element=err.Properties.VariableNames
    err=array2table(ones([1,length(vals.Properties.VariableNames)]).*0.1);
    err.Properties.VariableNames=vals.Properties.VariableNames;
    err.Spectrum={'Sigma'};
    signame={'Sigma'};
end
numelements=numel(vals.Properties.VariableNames);

%% DEBUG - plot corrolation plot and a basic plot
%{
figure()
scatter(vals.Ni,vals.Ta)
%}
h=figure();
%corrplot(vals)
corrplot(vals(:,varoffset:numelements-endoffset))
figname='CorrPlot';
saveas(h,[EDSfile(1:end-4),figname,'_base_',Baseelem],'tiffn')
saveas(h,[EDSfile(1:end-4),figname,'_base_',Baseelem],'fig')

%% rank sort
rank=table();
weightedvalue=rank;
for i=varoffset:numelements-endoffset
    
    elem=vals.Properties.VariableNames{i};
%     if strcmp(SEM,'Tescan')
%         elem=elemdat.Sym(strcmp(elemdat.Element,elem));
%     end
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

% assign solid fractions
vals.rankFs=(vals.ranksort-0.5)./size(vals,1);
vals.WIRSFs=(vals.WIRS-0.5)./size(vals,1);

%{%
%% find density
%first extract density of every element
elemdens=NaN(numelements,1);
for i=varoffset:numelements-endoffset
    elem=vals.Properties.VariableNames{i};
    if ~isempty(elemdat.Density(strcmp(elemdat.Sym,elem)))
        elemdens(i)=elemdat.Density(strcmp(elemdat.Sym,elem));
    elseif ~isempty(elemdat.Density(strcmp(elemdat.Element,elem)))
        elemdens(i)=elemdat.Density(strcmp(elemdat.Element,elem));
    end
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
saveas(h,[EDSfile(1:end-4),figname,'_base_',Baseelem],'tiffn')
saveas(h,[EDSfile(1:end-4),figname,'_base_',Baseelem],'fig')
%}
%% plot compositions as a function of Fs
h=figure();
scatter(vals.rankFs,table2array(vals(:,3:numelements-1)))
ylabel('Wt %')
xlabel('Fs')
legend(vals.Properties.VariableNames{3:numelements-1})
title('WIRS')
figname='RS Composition TOT';
saveas(h,[EDSfile(1:end-4),figname,'_base_',Baseelem],'tiffn')
saveas(h,[EDSfile(1:end-4),figname,'_base_',Baseelem],'fig')
%% plot each element individually on a subplot
h=figure();
for i=varoffset:numelements-endoffset
    subplot(3,4,i-varoffset+1)
    scatter(vals.WIRSFs,table2array(vals(:,i)),'r.')
    ylabel('Wt %')
    xlabel('Fs')
    title(vals.Properties.VariableNames{i})
end
sgtitle('WIRS')
figname='WIRS';
%saveas(h,[EDSfile(1:end-4),figname,'_base_',Baseelem,'.tif'],'tiffn')
print(h,[EDSfile(1:end-4),figname,'_V2_base_',Baseelem,'.tif'],'-dpng','-r1200');
saveas(h,[EDSfile(1:end-4),figname,'_base_',Baseelem,'.fig'],'fig')
%{
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
saveas(h,[EDSfile(1:end-4),figname,'_base_',Baseelem],'tiffn')
saveas(h,[EDSfile(1:end-4),figname,'_base_',Baseelem],'fig')
%}