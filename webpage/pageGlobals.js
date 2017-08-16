function add_new_color(){
  do{
    newCol='#'+Math.floor(Math.random()*16777215).toString(16);
  } while( colList.includes(newCol) );
  return newCol;//'#'+Math.floor(Math.random()*16777215).toString(16);
}

function sndDistChng(){
  d = document.getElementById("SelectADistrict");
  currentDistrict = d.options[d.selectedIndex].value;
}

function addAnotherDistrict(){
  others = [for (x of districts) if (x != 'Not Assigned') x];
  mycurrentDistrict = Math.max.apply(null, others) + 1;

  districts.push(mycurrentDistrict);
  allColors[mycurrentDistrict] = add_new_color();

  allDists = document.getElementById("SelectADistrict");
  oneToAdd = document.createElement("option");
  oneToAdd.text = mycurrentDistrict;
  allDists.options.add(oneToAdd, mycurrentDistrict);
}

function calculateAll(dist, funcName){
  var outputString = '';

  if(funcName == "compactness"){
    var totInterior = 0;

    map.data.forEach(function(feature){
      var val = feature.getProperty("District");

      if (val == dist){
        totInterior++;
      }
    });
    outputString += ("Compactness: " + totInterior.toString());
  }
  if(funcName == "contiguousness"){
    outputString += "Contiguousness: 1";
  }
  if(funcName == "area"){
    var totArea = 0;

    map.data.forEach(function(feature){
      var val = feature.getProperty("District");
      if (val == dist){
        totArea += feature.getProperty("ALAND10");
      }
    });
    outputString += ("Total Area: " + totArea.toFixed(3).toString());
  }
  return outputString;
}

function addToList(thingy){
  var values = document.getElementById("checkStats");
  if(thingy.checked){
    list_of_functions_to_compute.push(thingy.value);
  } else {
    list_of_functions_to_compute = [for (x of list_of_functions_to_compute) if (x != thingy.value) x];
  }
}

var colList = ["#000000"];
var jsonfilename="./NCPrecincts.json";
var connectionsfilename="../NorthCarolina/PRECINCTConnections.csv"
var centroid=[35.0, -79.9];
var map;
var currentDistrict="Not Assigned";
var districts=["Not Assigned",0];
var allColors={"Not Assigned":"#000000", 0:add_new_color()}

var yet_to_assign;
var list_of_functions_to_compute = [];


