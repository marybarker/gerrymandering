function add_new_color(){
  do{
    newCol = "rgb("+Math.floor(Math.random() * 255)+","+(Math.random() * 255)+","+(Math.random() * 255)+")";
  } while( colList.includes(newCol) );
  return newCol;
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
  document.getElementById("legendOfDists").innerHTML += '<li><span style="background:'+allColors[mycurrentDistrict]+'"></span>'+mycurrentDistrict.toString()+'</li>';
}

function calculateAll(dist, funcName){
  var outputString = '';

  if (funcName == "numInDist") {
    var tot = 0;
    map.data.forEach(function(feature){
      if (feature.getProperty("District") == dist) {
        tot++;
      }
    });
    outputString = ("Number of vtds in district: "+tot.toString()+"<br>");
  }
  /* COMPACTNESS */
  if(funcName == "compactness"){
    var totInterior = 0.0;
    var totExterior = 0.0;
    var vtds = [];
    map.data.forEach(function(feature){
      if (feature.getProperty("District") == dist){
        vtds.push(feature.getProperty("GEOID10"));
      }
    });
    
    for(var vtd in vtds){
      name1 = vtds[vtd];
      listA = [for (name2 of adjacencies[name1]) if (vtds.includes(name2)) name2];
      listB = [for (name2 of adjacencies[name1]) name2];
      if (listA.length == listB.length) {
        totInterior+=1.0;
      }else{
        totExterior+=1.0;
      }
    }
    if (totExterior < 1){
      totExterior = 1.0;
    }
    var cpctVal = totInterior / totExterior;
    outputString += ("Compactness: " + cpctVal.toFixed(3).toString()+"<br>");
  }
  /* CONTIGUOUSNESS */
  if(funcName == "contiguousness"){
    var contScore = 0;
    var vtds = [];

    map.data.forEach(function(feature){
      if (feature.getProperty("District") == dist){
        vtds.push(feature.getProperty("GEOID10"));
      }
    });
    while(vtds.length > 0){
      contScore++;
      var currentRegion = new Set();
      var addons = [ vtds[0] ];

      while( addons.length > 0) {
        for (i in addons){
          currentRegion.add(addons[i]);
        }
        var subsubedges = new Set([for (x of addons) for (y of adjacencies[x]) if (vtds.includes(y)) y]); 
        if (subsubedges.size > 0) {
          addons = Array.from(subsubedges);
          notToAdd = Array.from(currentRegion);
          addons = [for (x of addons) if (notToAdd.indexOf(x) < 0) x];
        }else{
          addons = [];
        }
      }
      vtds = [for (x of vtds) if (!(currentRegion.has(x))) x];
    }
    outputString += "Contiguousness: "+contScore+"<br>";
  }
  /* AREA */
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
  /* AFRICAN AMERICAN CONCENTRATION */
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
  /* HISPANIC CONCENTRATION */
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
  /* POPULATION */
  if(funcName == "population"){
    var totPop = 0;
    map.data.forEach(function(feature){
      if(feature.getProperty("District") == dist){
        totPop += parseFloat(blockstats[feature.getProperty("GEOID10")].population);
      }
    });
    outputString += ("Total Population: "+totPop.toString()+"<br>");
  }

  return outputString;
}

function updateCurrentStateInfo(){
  var numVTDS = calculateAll(currentDistrict, "numInDist");
  document.getElementById("currentDist").innerHTML="Current district: "+currentDistrict.toString();
  document.getElementById("unassigned").innerHTML="Total vtds not assigned: "+yet_to_assign.toString();

  var statsString='';
  for(var name of list_of_functions_to_compute){
    statsString += (calculateAll(currentDistrict, name)) + '<br>';
  }
  document.getElementById("stats").innerHTML = statsString;
}

function sndDistChng(){
  d = document.getElementById("SelectADistrict");
  currentDistrict = d.options[d.selectedIndex].value;
  updateCurrentStateInfo();
}

function addToList(thingy){
  if(thingy.checked){
    list_of_functions_to_compute.push(thingy.value);
  } else {
    list_of_functions_to_compute = [for (x of list_of_functions_to_compute) if (x != thingy.value) x];
  }
  updateCurrentStateInfo();
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

function loadAdjacencyFrame(data){
  var allRows = data.split("\n");
  var low, high, lenIdx;

  for(var singleRow=0; singleRow < allRows.length; singleRow++){
    var rowCells = allRows[singleRow].split(',');

    if (singleRow === 0) {
      lenIdx = rowCells.indexOf("length");
      low  = rowCells.indexOf("low");
      high = rowCells.indexOf("high");

    }else{
      var length = parseFloat(rowCells[lenIdx]);

      if (length > 0.0){
        if (rowCells[low] in adjacencies) {
          adjacencies[rowCells[low]].push(rowCells[high]);

        }else{
          adjacencies[rowCells[low]] = [rowCells[high]];

        }
        if (rowCells[high] in adjacencies) {
          adjacencies[rowCells[high]].push(rowCells[low]);

        }else{
          adjacencies[rowCells[high]] = [rowCells[low]];
        }
      }
    }
  }
}

function LoadStateAsDict(data){
  var allRows = data.split("\n");
  var vtdName, cdName, rowCells;
  numDistsForState = 0;
  
  for(var singleRow=0; singleRow < allRows.length; singleRow++){
    rowCells = allRows[singleRow].split(',');
    if (singleRow == 0) {
    }else{
      vtdName = rowCells[1];
      cdName = parseInt(rowCells[0]);
      currentState[vtdName] = cdName;
      if (cdName){
        numDistsForState = (numDistsForState > cdName) ? numDistsForState : cdName;
      }
    }
  }
}

function UseCurrentDistricting(){
  var createdDistricts = Math.max.apply(null, [for (x of districts) if (x != "Not Assigned") x]);

  for (var x = createdDistricts+1; x < (numDistsForState+1); x++){
    addAnotherDistrict();
    console.log(x);
  }
  
  map.data.forEach(function(feature){
    feature.setProperty("District", currentState[feature.getProperty("GEOID10")]);
  });
  yet_to_assign = 0;
  updateCurrentStateInfo();
}

/* Don't need generalizing */
var map;
var yet_to_assign;
var list_of_functions_to_compute = [];

/* depends on which state you want */
var jsonfilename="./NC/NCPrecincts.json";
var connectionsfilename="./NC/PRECINCTconnections.csv"
var statsfilename="./NC/vtdstats.csv"
var cdsfilename="./NC/CurrentState.csv"
var centroid=[35.0, -79.9];
var indexingCol = "GEOID10";
var listOfCols = ["ID", "ALAND", "aframcon", "hispcon", "population"];
var blockstats = {};
var adjacencies = {};
var currentState={};
var numDistsForState;

/* depend on the current coding setup */
var colList = ["#000000"];
var currentDistrict="Not Assigned";
var districts=["Not Assigned",0];
var allColors={"Not Assigned":"#000000", 0:add_new_color()};


$.ajax({
  url:statsfilename,
  dataType:'text',
}).done(loadCSVtoArrays)

$.ajax({
  url:connectionsfilename,
  dataType:'text',
}).done(loadAdjacencyFrame) 

$.ajax({
  url:cdsfilename,
  dataType:'text',
}).done(LoadStateAsDict)


