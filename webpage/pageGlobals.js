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
    outputString += ("Total Area: " + totArea.toFixed(3).toString()+'<br>');
  }
  if(funcName == "aframcon"){
    var totConc = 0.0;
    var toDivide = 0;
    var idName, addValue,divValue;

    map.data.forEach(function(feature){
      var val = feature.getProperty("District");
      if (val == dist){
        idName = feature.getProperty("GEOID10");

        divValue = parseFloat(blockstats[idName].population);
        addValue = parseFloat(blockstats[idName].aframcon);

        toDivide += 1.0 * divValue;
        totConc += (divValue * addValue);
      }
    });
    totConc /= toDivide;
    outputString += ("African American Concentration: "+totConc.toFixed(3).toString()+'<br>');
  }
  if(funcName == "hispcon"){
    var totConc = 0.0;
    var toDivide = 0;
    var idName, addValue,divValue;

    map.data.forEach(function(feature){
      var val = feature.getProperty("District");
      if (val == dist){
        idName = feature.getProperty("GEOID10");

        divValue = parseFloat(blockstats[idName].population);
        addValue = parseFloat(blockstats[idName].hispcon);

        toDivide += 1.0 * divValue;
        totConc += (divValue * addValue);
      }
    });
    totConc /= toDivide;
    outputString += ("Hispanic Concentration: "+totConc.toFixed(3).toString()+'<br>');
  }
  return outputString;
}

function addToList(thingy){
  if(thingy.checked){
    list_of_functions_to_compute.push(thingy.value);
  } else {
    list_of_functions_to_compute = [for (x of list_of_functions_to_compute) if (x != thingy.value) x];
  }
}

function zip(a, b){
  var arr = {};
  for(var key in a){
    arr[a[key]] = b[key];
  } 
  return arr;
}

function loadCSVtoArrays(data){
  var allRows = data.split("\n");
  var idIdx;
  var toReturn;
  var otherList = [];
  
  for(var singleRow=0; singleRow < allRows.length; singleRow++){
    var rowCells = allRows[singleRow].split(',');
    
    if (singleRow === 0){
      idIdx = rowCells.indexOf(indexingCol);
      otherList = [for (x of listOfCols) rowCells.indexOf(x)];
    }else{
      blockstats[rowCells[idIdx]] = zip(listOfCols, [for (x of otherList) rowCells[x]]);
    }
  }
}

/* Don't need generalizing */
var map;
var yet_to_assign;
var list_of_functions_to_compute = [];

/* depends on which state you want */
var jsonfilename="./NC/NCPrecincts.json";
var connectionsfilename="./NC/PRECINCTconnections.csv"
var statsfilename="./NC/vtdstats.csv"
var centroid=[35.0, -79.9];
var indexingCol = "GEOID10";
var listOfCols = ["ID", "ALAND", "aframcon", "hispcon", "population"];
var blockstats = {};

/* depend on the current coding setup */
var colList = ["#000000"];
var currentDistrict="Not Assigned";
var districts=["Not Assigned",0];
var allColors={"Not Assigned":"#000000", 0:add_new_color()};

$.ajax({
  url:statsfilename,
  dataType:'text',
}).done(loadCSVtoArrays)

//console.log(blockstats);


