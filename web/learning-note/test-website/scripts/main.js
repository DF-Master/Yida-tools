let myHeading = document.querySelector('h1');
myHeading.onclick = function () {
	alert('别戳我，我怕疼。');
};

let myImage = document.querySelector('img');

myImage.onclick = function () {
	let mySrc = myImage.getAttribute('src');
	if (mySrc === 'images/q1.png') {
		myImage.setAttribute('src', 'images/q2.png');
	} else {
		myImage.setAttribute('src', 'images/q1.png');
	}
};

let myButton = document.querySelector('button');

myButton.onclick = function () {
	setUserName();
};

function setUserName() {
	let myName = prompt('请输入你的名字。');
	if (!myName) {
		setUserName();
	} else {
		localStorage.setItem('name', myName);
		myHeading.textContent = 'Hello ' + myName;
	}
}

if (!localStorage.getItem('name')) {
	setUserName();
} else {
	let storedName = localStorage.getItem('name');
	myHeading.textContent = 'Hello ' + storedName;
}
