function machines=getMachines(NumOfMachines)
%% Usage: this function will return the machines string based on the total number of machines specified in the input
%% Input: the number of machines
%% Output: machines string
%%

switch NumOfMachines
           case 1
                  machines={'n117'};
           case 2
                  machines={'n117' 'n118'};
           case 3
                  machines={'n117' 'n118' 'n119'};
           case 4
                  machines={'n117' 'n118' 'n119' 'n120'};
           case 5
                  machines={'n117' 'n118' 'n119' 'n120' 'n121'};
           case 6
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122'};
           case 7
                  machines={'n117' 'n118' 'n119' 'n121' 'n122' 'n123' 'n124'};
           case 8
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124'};
           case 9
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125'};
           case 10
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126'};
          case 11
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127'};
          case 12
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128'};
          case 13
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129'};
          case 14
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130'};
          case 15
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131'};
          case 16
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131' 'n132'};
          case 17 
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131' 'n132' 'n101'};
          case 18
		  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131' 'n132' 'n101' 'n102'};
          case 19
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131' 'n132' 'n101' 'n102' 'n103'};
          case 32 
		  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131' 'n132' 'n101' 'n102' 'n103' 'n104' 'n105' 'n106' 'n107' 'n108' 'n109' 'n110' 'n111' 'n112' 'n113' 'n114' 'n115' 'n116' };
	  case 31
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131' 'n132' 'n101' 'n102' 'n103' 'n104' 'n105' 'n106' 'n107' 'n108' 'n109' 'n110' 'n111' 'n112' 'n113' 'n114' 'n115' };
	  case 30
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131' 'n132' 'n101' 'n102' 'n103' 'n104' 'n105' 'n106' 'n107' 'n108' 'n109' 'n110' 'n111' 'n112' 'n113' 'n114' };
          case 29
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131' 'n132' 'n101' 'n102' 'n103' 'n104' 'n105' 'n106' 'n107' 'n108' 'n109' 'n110' 'n111' 'n112' 'n113' };	  
 	  case 28
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131' 'n132' 'n101' 'n102' 'n103' 'n104' 'n105' 'n106' 'n107' 'n108' 'n109' 'n110' 'n111' 'n112' };
          case 27
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131' 'n132' 'n101' 'n102' 'n103' 'n104' 'n105' 'n106' 'n107' 'n108' 'n109' 'n110' 'n111' };
          case 26
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131' 'n132' 'n101' 'n102' 'n103' 'n104' 'n105' 'n106' 'n107' 'n108' 'n109' 'n110'};
          case 25
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131' 'n132' 'n101' 'n102' 'n103' 'n104' 'n105' 'n106' 'n107' 'n108' 'n109' };
          case 24
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131' 'n132' 'n101' 'n102' 'n103' 'n104' 'n105' 'n106' 'n107' 'n108' };
          case 23
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131' 'n132' 'n101' 'n102' 'n103' 'n104' 'n105' 'n106' 'n107' };
	  case 22
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131' 'n132' 'n101' 'n102' 'n103' 'n104' 'n105' 'n106' };
          case 21
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131' 'n132' 'n101' 'n102' 'n103' 'n104' 'n105' };
          case 20
                  machines={'n117' 'n118' 'n119' 'n120' 'n121' 'n122' 'n123' 'n124' 'n125' 'n126' 'n127' 'n128' 'n129' 'n130' 'n131' 'n132' 'n101' 'n102' 'n103' 'n104' };
end
