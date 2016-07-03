for pubCount=1:8;
    LIST=[0.2 0.5 0.8 1 1.2 1.4 5 10];
    MODIFIER=LIST(pubCount);
    save('temp.mat','MODIFIER','pubCount');
    
publish('gg_EnergySinkSimulator.m');
publish('gg_EnergySinkGrapher.m');
publish('gg_EnergySinkTyper.m');
system(['ren html ',simname]);
system(['move ' ,simname,'*.* ', simname, '\']);

load('temp.mat');
end;
