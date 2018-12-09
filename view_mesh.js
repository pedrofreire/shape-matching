
let scenewrapper;

function start() {
  scenewrapper = new SceneWrapper();
  const loader = new THREE.OBJLoader();
  loader.load('./datasets/non-rigid-world/cat0.obj', (obj) => { 
  // loader.load('./datasets/simple/reg_tetra.obj', (obj) => { 
    scenewrapper.scene.add(obj); 
  });
}

function display_sample() {
  const loader = new THREE.FileLoader();
  loader.load('./sample4.vert', (data) => {
    const lines = data.split('\n');
    console.log(lines);
    for(let line of lines) {
      const nums = parse_nums(line);
      add_sphere(nums);
    }
  });
}

function parse_nums(line) {
  console.log(line);
  const str_nums = line.split(' ');
  console.log(str_nums);
  return str_nums.map(x => parseFloat(x));
}

function add_sphere(nums) {
  const [x, y, z] = nums;
  console.log(nums);
  const radius = 1;
  const resolution = 30;
  const geometry = new THREE.SphereGeometry(radius);
  const material = new THREE.MeshBasicMaterial({color : 0xEEEE00});

  const sphere = new THREE.Mesh(geometry, material);
  sphere.position.set(x, y, z);

  scenewrapper.scene.add(sphere);
}

window.onload = function () {
  start();
  display_sample();
}
