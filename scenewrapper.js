
class SceneWrapper {
	constructor() {
		this.initCamera();
		this.initRenderer();
		this.initScene();
		this.animations = [];
		this.animate(0);
	}

	initCamera() {
		const angle = 60;
		const ratio = window.innerWidth/window.innerHeight;
		const nearClippingDist = 1;
		const farClippingDist = 2000;
		this.camera = new THREE.PerspectiveCamera(angle, ratio, nearClippingDist, farClippingDist);
		this.camera.position.x = 100;
		this.camera.position.y = -200;
		this.camera.position.z = 150;
		this.camera.lookAt(new THREE.Vector3(0,0,0));
	}

	initRenderer() {
		this.renderer = new THREE.WebGLRenderer();
		this.renderer.setSize(window.innerWidth, window.innerHeight);
		document.body.appendChild(this.renderer.domElement);
	}

	initScene() {
		this.scene = new THREE.Scene();
		this.orbitControls = new THREE.OrbitControls(this.camera);
		this.initLight();
	}

	render() {
		this.renderer.render(this.scene, this.camera);
	}

	initLight() {
        var spotLight = new THREE.SpotLight(0xffffff,0.9);
        var lightPosition = new THREE.Vector3(0, 100, 200);
        spotLight.position.copy(lightPosition);
        spotLight.castShadow = true;
        this.scene.add(spotLight);

        spotLight = new THREE.SpotLight(0xffffff,0.9);
        lightPosition = new THREE.Vector3(-30, 100, 30);
        spotLight.position.copy(lightPosition);
        spotLight.castShadow = true;
        this.scene.add(spotLight);

        spotLight = new THREE.SpotLight(0xffffff,0.9);
        lightPosition = new THREE.Vector3(200, 200, 200);
        spotLight.position.copy(lightPosition);
        spotLight.castShadow = true;
        this.scene.add(spotLight);

	    const ambient = new THREE.AmbientLight( 0xffffff, 0.5 );
	    this.scene.add(ambient);
	}
	
	getAllObjects(obj) {
		if(obj == undefined)
			obj = this.scene;

		var objects = [obj];
		for(let child of obj.children)
			objects.push(...this.getAllObjects(child));

		return objects;
	}

	animate(timeEllapsed) {
		requestAnimationFrame(this.animate.bind(this));
		// Input.handle();
		for(let animation of this.animations)
			animation(timeEllapsed);
		this.render();
	}
}

