// game.js
const canvas = document.getElementById('game');
const ctx = canvas.getContext('2d');

let snake = [
	{ x: 150, y: 150 },
	{ x: 140, y: 150 },
	{ x: 130, y: 150 },
	{ x: 120, y: 150 },
	{ x: 110, y: 150 },
];
let dx = 10;
let dy = 0;
let foodX;
let foodY;
let score = 0;

function draw() {
	ctx.clearRect(0, 0, canvas.width, canvas.height);

	// Draw the snake
	snake.forEach((segment) => {
		ctx.fillStyle = 'lightgreen';
		ctx.strokeStyle = 'darkgreen';
		ctx.fillRect(segment.x, segment.y, 10, 10);
		ctx.strokeRect(segment.x, segment.y, 10, 10);
	});

	// Draw the food
	ctx.fillStyle = 'red';
	ctx.strokeStyle = 'darkred';
	ctx.fillRect(foodX, foodY, 10, 10);
	ctx.strokeRect(foodX, foodY, 10, 10);

	// Move the snake
	const head = { x: snake[0].x + dx, y: snake[0].y + dy };
	snake.unshift(head);
	const ateFood = snake[0].x === foodX && snake[0].y === foodY;
	if (ateFood) {
		score += 10;
		document.getElementById('score').innerHTML = score;
		generateFood();
	} else {
		snake.pop();
	}

	// Check for game over
	checkGameOver();
}

function checkGameOver() {
	// Check if the snake has collided with the wall
	if (
		snake[0].x < 0 ||
		snake[0].x > canvas.width - 10 ||
		snake[0].y < 0 ||
		snake[0].y > canvas.height - 10
	) {
		gameOver();
		return;
	}

	// Check if the snake has collided with itself
	for (let i = 1; i < snake.length; i++) {
		if (snake[0].x === snake[i].x && snake[0].y === snake[i].y) {
			gameOver();
			return;
		}
	}
}

function gameOver() {
	// Stop the game loop
	clearInterval(gameTimer);

	// Show the game over message
	ctx.font = '48px sans-serif';
	ctx.fillStyle = 'red';
	ctx.textAlign = 'center';
	ctx.fillText('Game Over', canvas.width / 2, canvas.height / 2);
}

function generateFood() {
	foodX = Math.floor(Math.random() * ((canvas.width - 10) / 10)) * 10;
	foodY = Math.floor(Math.random() * ((canvas.height - 10) / 10)) * 10;

	// Check if the food is generated on top of the snake
	snake.forEach((segment) => {
		if (segment.x === foodX && segment.y === foodY) generateFood();
	});
}

function changeDirection(event) {
	const LEFT_KEY = 37;
	const RIGHT_KEY = 39;
	const UP_KEY = 38;
	const DOWN_KEY = 40;

	const keyPressed = event.keyCode;
	const goingUp = dy === -10;
	const goingDown = dy === 10;
	const goingRight = dx === 10;
	const goingLeft = dx === -10;

	if (keyPressed === LEFT_KEY && !goingRight) {
		dx = -10;
		dy = 0;
	}
	if (keyPressed === UP_KEY && !goingDown) {
		dx = 0;
		dy = -10;
	}
	if (keyPressed === RIGHT_KEY && !goingLeft) {
		dx = 10;
		dy = 0;
	}
	if (keyPressed === DOWN_KEY && !goingUp) {
		dx = 0;
		dy = 10;
	}
}

const startButton = document.getElementById('start-button');
const pauseButton = document.getElementById('pause-button');

startButton.addEventListener('click', startGame);
pauseButton.addEventListener('click', pauseGame);

function startGame() {
	document.addEventListener('keydown', changeDirection);
	gameTimer = setInterval(draw, 100);
	// code to start the game
}

function pauseGame() {
	clearInterval(gameTimer);
	// code to pause the game
}

generateFood();
draw();
