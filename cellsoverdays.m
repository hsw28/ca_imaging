function [repeatnumsplacecells repeatnumsallcells fwd_field bwd_field] = cellsoverdays(alignmentdata, fieldcentersSTRUCTURE)
%computes number of cells that repeat over days and which days
%also gives field centers for repeated cells
%structure should be made from CIATAH

numses = NaN(length(alignmentdata),1);

daynames = fieldnames(fieldcentersSTRUCTURE);
daynum = length(daynames);

fwd_field = NaN(length(alignmentdata),6);
bwd_field = NaN(length(alignmentdata),6);

repeatnumsplacecellsfwd = [];
repeatnumsplacecellsbwd = [];
repeatnumsplacecellsall = [];
figure

for k=1:size(alignmentdata, 1)
  numses(k) = length(find(alignmentdata(k,:)>0)); %find number of sessions
  if length(find(alignmentdata(k,:)>0))>0
    for z=1:length(alignmentdata(k,:)) %z is the day
      if alignmentdata(k,z)>0
        name = char(daynames(z));
        currentday = fieldcentersSTRUCTURE.(name);
        fwd_field(k,z) = currentday(alignmentdata(k,z),1);
        bwd_field(k,z) = currentday(alignmentdata(k,z),2);
      end
    end
%{
    subplot(3,2,length(find(alignmentdata(k,:)>0))-1)

    if length(find(isnan(fwd_field(k,:))==0)) ==length(find(alignmentdata(k,:)>0))
      plot(fwd_field(k,:))
        hold on
    end
    if length(find(isnan(bwd_field(k,:))==0))==length(find(alignmentdata(k,:)>0))
      plot(bwd_field(k,:))
        hold on
    end
%}


  end
  all = nansum([fwd_field(k,:);bwd_field(k,:)],1);
  repeatnumsplacecellsfwd(end+1) = length(find(isnan(fwd_field(k,:))==0));
  repeatnumsplacecellsbwd(end+1) = length(find(isnan(bwd_field(k,:))==0));
  repeatnumsplacecellsall(end+1) = length(find(all>0));
end

repeatnumsplacecells = [repeatnumsplacecellsfwd;repeatnumsplacecellsbwd; repeatnumsplacecellsall]';
repeatnumsplacecells = [histcounts(repeatnumsplacecells(:,1), [-05,.5,1.5,2.5,3.5,4.5,5.5,6.5]); histcounts(repeatnumsplacecells(:,2), [-05,.5,1.5,2.5,3.5,4.5,5.5,6.5]); histcounts(repeatnumsplacecells(:,3),[-05,.5,1.5,2.5,3.5,4.5,5.5,6.5])];
repeatnumsplacecells = (repeatnumsplacecells./length(numses)).*100;
repeatnums= numses;
f = repeatnums;
per = [];
per(end+1) = length(find(f==1))./length(f);
per(end+1) = length(find(f==2))./length(f);
per(end+1) = length(find(f==3))./length(f);
per(end+1) = length(find(f==4))./length(f);
per(end+1) = length(find(f==5))./length(f);
per(end+1) = length(find(f==6))./length(f);
per;
