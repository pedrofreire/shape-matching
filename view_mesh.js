

function start() {
  const scenewrapper = new SceneWrapper();
  const loader = new THREE.OBJLoader();
  loader.load('./datasets/non-rigid-world/cat1.obj', (obj) => { 
    scenewrapper.scene.add(obj); 
  });
}

window.onload = function () {
  start();
}
