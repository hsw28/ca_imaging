function f = local_linearregression(localI, pos, spiketimes, fieldcenters);

  %speed in place cell bin
  %time in bin
  %fields in bin
  %distance from reward
  %spike rate


  if isstruct(fieldcenters)==1
    fieldcenter_names = fieldnames(fieldcenters);
    localI_names = fieldnames(localI);
    pos_names = fieldnames(pos);
    spiketimes_names =  fieldnames(spiketimes);
    daynum = length(fieldcenter_names);
  else
    daynum =1;
  end

  x= [];
  y = [];
  for z=1:(daynum)
    if isstruct(fieldcenters)==1
      name = char(fieldcenter_names(z))
      currentfieldcenter = fieldcenters.(name);

      name = char(localI_names(z))
      currentlocalI = localI.(name);

      name = char(pos_names(z))
      currentpos = pos.(name);

      name = char(spiketimes_names(z))
      currentspiketimes = spiketimes.(name);

      velcorr= [];
      timeinfield= [];
      distfromrew= [];
      eventrate_all= [];
      eventrate_moving= [];
      maxrate= [];

    else
      currentfieldcenter = fieldcenters;
      currentlocalI = localI;
      currentpos = pos;
      currentspiketimes = spiketimes;
    end

  temp = [currentlocalI(:,1); currentlocalI(:,3)].*100;

  y = [y; temp];
  currenty = temp;

  velcorr = localIvsvel(currentlocalI, currentpos);
  velcorr = velcorr(:,1);


  timeinfield = localIvstime(currentlocalI, currentpos);
  timeinfield = timeinfield(:,1);

  distfromrew = [264-currentlocalI(:,2).*4; currentlocalI(:,4).*4];
  length(distfromrew)

  want = length(currenty)./2;
  empty = NaN(1,want);
  temp1 = find(currentfieldcenter(:,1)>0);
  temp2 = find(currentfieldcenter(:,2)>0);
  eventrate = ca_firingrate(currentspiketimes, currentpos);
  eventrate1 = eventrate(temp1);
  eventrate2 = eventrate(temp2);
  r1 = empty;
  r1(1:length(eventrate1))= eventrate1;
  r2 = empty;
  r2(1:length(eventrate2))= eventrate2;
  eventrate_all = ([r1, r2])';


  maxrate = ca_maxrate(currentspiketimes, currentpos, 4);
  meanrate1 = maxrate(:,3);
  meanrate2 = maxrate(:,4);
  meanrate1 = meanrate1(temp1);
  meanrate2 = meanrate2(temp2);
  m3 = empty;
  m3(1:length(meanrate1))= meanrate1;
  m4 = empty;
  m4(1:length(meanrate2))= meanrate2;
  eventrate_moving = ([m3, m4])';

  maxrate1 = maxrate(:,1);
  maxrate2 = maxrate(:,2);
  maxrate1 = maxrate1(temp1);
  maxrate2 = maxrate2(temp2);
  m1 = empty;
  m1(1:length(maxrate1))= maxrate1;
  m2 = empty;
  m2(1:length(maxrate2))= maxrate2;
  maxrate = ([m1, m2])';


  size(velcorr)
  size(timeinfield)
  size(distfromrew)

  size(eventrate_all)
  size(eventrate_moving)
  size(maxrate)
  temp = [velcorr, timeinfield, distfromrew, eventrate_all, eventrate_moving, maxrate];
  x = [x; temp];


[beta,Sigma,E,CovB,logL] = mvregress(x, y);
mdl = fitlm(x, y);

end

mdl = fitlm(x, y)
