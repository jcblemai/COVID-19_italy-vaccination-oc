Spro(9).DEN_PROV='ViboValentia';
Spro(10).DEN_PROV='VerbanoCusioOssola';
Spro(11).DEN_PROV='MonzaBrianza';
Spro(13).DEN_PROV='BarlettaAndriaTrani';
Spro(14).DEN_PROV='Torino';
Spro(23).DEN_PROV='Genova';
Spro(24).DEN_PROV='LaSpezia';
Spro(28).DEN_PROV='Milano';
Spro(40).DEN_PROV='Venezia';
Spro(48).DEN_PROV='ReggioEmilia';
Spro(50).DEN_PROV='Bologna';
Spro(53).DEN_PROV='ForliCesena';
Spro(54).DEN_PROV='PesaroeUrbino';
Spro(57).DEN_PROV='AscoliPiceno';
Spro(58).DEN_PROV='MassaCarrara';
Spro(61).DEN_PROV='Firenze';
Spro(71).DEN_PROV='Roma';
Spro(76).DEN_PROV='Napoli'; % ex medio campidarno
Spro(79).DEN_PROV='LAquila';
Spro(85).DEN_PROV='Bari';
Spro(93).DEN_PROV='ReggioCalabria';
Spro(95).DEN_PROV='Palermo';
Spro(96).DEN_PROV='Messina';
Spro(100).DEN_PROV='Catania';
Spro(105).DEN_PROV='Cagliari';
Spro(107).DEN_PROV='SudSardegna';

% find key between provinces in NodeCell and shapefile
shpKey=zeros(length(NodeCell),1);
for i=1:length(NodeCell)
    for j=1:length(Spro)
        if strcmp(NodeCell{i,1},Spro(j).DEN_PROV)
            shpKey(i)=j;
            break
        end
    end
end
