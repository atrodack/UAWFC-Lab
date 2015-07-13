function [Ppos_Full] = setOverlapActs(pistonlist, slaveActs)



for n = 1:length(slaveActs)
    switch slaveActs(n)
        case 572
            
            NearbyActs = [pistonlist(540) pistonlist(539) pistonlist(571) pistonlist(603)];
            pistonlist(slaveActs(n),1) = mean(NearbyActs);
            
            
        case 635
            
            NearbyActs = [pistonlist(603) pistonlist(602) pistonlist(634)];
            pistonlist(slaveActs(n),1) = mean(NearbyActs);
            
        case 666
            
            NearbyActs = [pistonlist(633) pistonlist(634) pistonlist(665) pistonlist(697) pistonlist(635)];
            pistonlist(slaveActs(n),1) = mean(NearbyActs);
            
        case 667
            
            NearbyActs = [pistonlist(635) pistonlist(634) pistonlist(666) pistonlist(635)];
            pistonlist(slaveActs(n),1) = mean(NearbyActs);
                        
        case 698
                        
            NearbyActs = [pistonlist(697) pistonlist(665) pistonlist(667) pistonlist(666)];
            pistonlist(slaveActs(n),1) = mean(NearbyActs);            
                        
        case 729
                        
            NearbyActs = [pistonlist(697) pistonlist(696) pistonlist(698) pistonlist(728) pistonlist(760)];
            pistonlist(slaveActs(n),1) = mean(NearbyActs);
                        
        case 730          
            
            NearbyActs = [pistonlist(697) pistonlist(729) pistonlist(698)];
            pistonlist(slaveActs(n),1) = mean(NearbyActs);
                        
        case 761
                        
            NearbyActs = [pistonlist(730) pistonlist(729) pistonlist(728) pistonlist(760)];
            pistonlist(slaveActs(n),1) = mean(NearbyActs);           
            
        case 792
            
            NearbyActs = [pistonlist(759) pistonlist(761) pistonlist(791) pistonlist(760)];                       
            pistonlist(slaveActs(n),1) = mean(NearbyActs);
                        
        case 823
                        
            NearbyActs = [pistonlist(790) pistonlist(791) pistonlist(792) pistonlist(822)];
            pistonlist(slaveActs(n),1) = mean(NearbyActs);
            
        otherwise
            error('You have Changed the Overlap Region! This must be re-calibrated!');
    end
end

Ppos_Full = pistonlist;











end